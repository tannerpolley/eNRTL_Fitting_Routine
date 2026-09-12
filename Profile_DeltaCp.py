import os

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyomo.environ as pyo
from scipy.linalg import eigh
from idaes.core.solvers import get_solver
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models_extra.column_models.properties import ModularPropertiesInherentReactionsInitializer

from Fitting_Routine import column_names, get_mole_fraction, param_dic, species_dic
from Load_Datasets import add_ABS_dataset, load_datasets, heat_of_absorption_expression
from Parameter_Setup import get_estimated_params, load_fitted_params, setup_param_scaling
from Uncertainty_Analysis import get_reduced_hessian
from eNRTL_property_setup import get_prop_dict


DATA_FILE = os.path.join("data", "Plots", "DeltaCp_Profile.csv")
PLOT_FILE = os.path.join("data", "Plots", "DeltaCp_Profile.png")
DATASET_DIR = os.path.join("data", "data_sets_to_load")
FIT_PARAMETER_FILE = os.path.join("data", "Parameters", "Parameters_fit.csv")
CANDIDATE_PARAMETER_FILE = os.path.join(
    "data", "Parameters", "Parameters_historical_coordinate.csv"
)
CANDIDATE_CORRELATION_FILE = os.path.join(
    "data", "Parameters", "Parameter_Correlation_historical_coordinate.csv"
)
CANDIDATE_SUMMARY_FILE = os.path.join(
    "data", "Plots", "Historical_Coordinate_Summary.csv"
)
REGULARIZATION = 0.01

SOLVER_OPTIONS = {
    "linear_solver": "ma57",
    "OF_ma57_automatic_scaling": "yes",
    "max_iter": 300,
    "tol": 1e-8,
    "constr_viol_tol": 1e-8,
    "halt_on_ampl_error": "no",
}


def metrics(values):
    values = np.asarray(values, dtype=float)
    return {
        "mae": float(np.mean(np.abs(values))),
        "rmse": float(np.sqrt(np.mean(values**2))),
        "bias": float(np.mean(values)),
    }


def heat_metrics(model, block, data, block_start=0):
    residuals = []
    start = block_start
    group_columns = (
        (["source"] if "source" in data.columns else [])
        + [column_names["temperature"], "experiment"]
    )
    for _, group in data.groupby(group_columns, sort=False):
        old = block[start]
        for offset, (_, row) in enumerate(group.iterrows(), start=1):
            current = block[start + offset]
            predicted = pyo.value(heat_of_absorption_expression(model, old, current)) * 1e-3
            residuals.append(predicted - row[column_names["heat_of_absorption"]])
            old = current
        start += len(group) + 1
    return metrics(residuals)


def build_model():
    model = pyo.ConcreteModel()
    model.params = GenericParameterBlock(**get_prop_dict(species_dic["components"]))
    setup_param_scaling(model)
    objective = 0
    objective, _, block_names = load_datasets(
        model,
        objective,
        DATASET_DIR,
        species_dic,
        get_mole_fraction,
        column_names,
        exclude_list=["Xu", "Bottinger"],
    )

    iscale.calculate_scaling_factors(model)
    initializer = ModularPropertiesInherentReactionsInitializer(
        solver="ipopt", solver_options=SOLVER_OPTIONS, output_level=idaeslog.WARNING
    )
    for block_name in block_names:
        initializer.initialize(getattr(model, block_name))

    unfit = get_estimated_params(model, param_dic)
    variables = list(unfit["Object"])
    reference_values = {var.name: pyo.value(var) for var in variables}
    fitted = pd.read_csv(FIT_PARAMETER_FILE)
    load_fitted_params(model, fitted)

    validation = pd.read_csv(os.path.join(DATASET_DIR, "kim_2014_dHabs.csv"))
    model.validation_states = model.params.build_state_block(
        range(len(validation) + 3), defined_state=True
    )
    add_ABS_dataset(
        model,
        model.validation_states,
        validation,
        column_names,
        species_dic,
        get_mole_fraction,
        0,
    )
    iscale.calculate_scaling_factors(model)
    initializer.initialize(model.validation_states)

    objective += sum(
        REGULARIZATION
        * 0.5
        * ((var - reference_values[var.name]) * iscale.get_scaling_factor(var)) ** 2
        for var in variables
    )
    model.obj = pyo.Objective(expr=objective)
    scaled = pyo.TransformationFactory("core.scale_model").create_using(
        model, rename=False
    )
    scaled.ipopt_zL_out = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    scaled.ipopt_zU_out = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    scaled_variables = [scaled.find_component(var.name) for var in variables]
    for var in scaled_variables:
        var.unfix()

    fit_blocks = []
    for filename in sorted(os.listdir(DATASET_DIR)):
        name, year, dataset_type = filename.split("_")
        dataset_type = dataset_type.split(".")[0]
        if name in {"Xu", "Bottinger"}:
            continue
        data = pd.read_csv(os.path.join(DATASET_DIR, filename))
        if dataset_type == "dHabs":
            data = data[data["fit"].astype(bool)]
            if data.empty:
                continue
        fit_blocks.append(
            (dataset_type, data, getattr(model, f"{name}_{year}_{dataset_type}"))
        )

    return model, scaled, variables, scaled_variables, fitted, fit_blocks, validation


def equilibrium_metrics(fit_blocks):
    """Pool observation residuals, retaining each source and explicit units."""
    records = []
    for kind, data, block in sorted(fit_blocks, key=lambda item: item[2].local_name):
        source = block.local_name
        for i, (_, row) in enumerate(data.iterrows()):
            if kind == "VLE" and 0.1 < row[column_names["loading"]] < 0.6:
                observed = float(row[column_names["CO2_pressure"]])
                predicted = pyo.value(block[i].fug_phase_comp["Liq", "CO2"]) / 1e3
                records.append(dict(kind="vle", source=source, row=i, species="CO2",
                                    observed=observed, predicted=predicted, unit="kPa"))
            elif kind == "ChEq":
                excluded = {column_names["amine_concentration"], column_names["temperature"],
                            column_names["loading"], "CO2", "CO3^2-"}
                for species in sorted(set(data.columns) - excluded):
                    records.append(dict(kind="speciation", source=source, row=i,
                                        species=species, observed=float(row[species]),
                                        predicted=pyo.value(block[i].mole_frac_phase_comp_true["Liq", species]),
                                        unit="mol/mol true species"))
    frame = pd.DataFrame(records)
    result = {}
    for kind, group in frame.groupby("kind", sort=True):
        def summarize(rows):
            residual = rows.predicted.to_numpy() - rows.observed.to_numpy()
            summary = dict(n=len(rows), unit=rows.unit.iloc[0], **metrics(residual))
            if kind == "vle":
                fractional = residual / rows.observed.to_numpy()
                summary.update(log_rmse=metrics(np.log(rows.predicted / rows.observed))["rmse"],
                               pressure_mape_fraction=float(np.mean(np.abs(fractional))),
                               pressure_mape_percent=float(100 * np.mean(np.abs(fractional))))
            return summary
        result[kind] = summarize(group)
        result[kind]["sources"] = {
            source: summarize(rows) for source, rows in group.groupby("source", sort=True)
        }
    return result, frame

def evaluate(model, fit_blocks, validation):
    result, _ = equilibrium_metrics(fit_blocks)
    for dataset_type, data, block in fit_blocks:
        if dataset_type == "dHabs":
            result["fit_heat"] = heat_metrics(model, block, data)
    result["external"] = heat_metrics(model, model.validation_states, validation)
    validation_groups = list(
        validation.groupby(
            (["source"] if "source" in validation.columns else [])
            + [column_names["temperature"], "experiment"],
            sort=False,
        )
    )
    validation_starts = {}
    start = 0
    for key, group in validation_groups:
        validation_starts[key] = start
        start += len(group) + 1
    for temperature, group in validation.groupby("temperature"):
        temperature_start = min(
            start
            for key, start in validation_starts.items()
            if temperature == (key[1] if "source" in validation.columns else key[0])
        )
        result[f"external_{int(temperature)}"] = heat_metrics(
            model,
            model.validation_states,
            group,
            block_start=temperature_start,
        )
    return result


def set_parameter_start(scaled_variables, fitted, scales):
    fitted_values = dict(zip(fitted["Object_Name"], fitted["Value"]))
    for var in scaled_variables:
        var.unfix()
        var.set_value(fitted_values[var.name] * scales[var.name])


def run_profile():
    model, scaled, variables, scaled_variables, fitted, fit_blocks, validation = build_model()
    solver = get_solver("ipopt", options=SOLVER_OPTIONS)
    fitted_values = dict(zip(fitted["Object_Name"], fitted["Value"]))
    specifications = {
        "bicarbonate": (
            "params.reaction_MEA_bicarbonate_formation_combo.dcp_rxn",
            float(fitted_values["params.reaction_MEA_bicarbonate_formation_combo.dcp_rxn"]),
            float(fitted.set_index("Object_Name").loc["params.reaction_MEA_bicarbonate_formation_combo.dcp_rxn", "Curvature_Scale"]),
        ),
        "carbamate": (
            "params.reaction_MEA_carbamate_formation_combo.dcp_rxn",
            float(fitted_values["params.reaction_MEA_carbamate_formation_combo.dcp_rxn"]),
            float(fitted.set_index("Object_Name").loc["params.reaction_MEA_carbamate_formation_combo.dcp_rxn", "Curvature_Scale"]),
        ),
    }
    scales = {
        var.name: iscale.get_scaling_factor(var, default=1) for var in variables
    }
    records = []
    baseline_objective = None

    for reaction, (target_name, center, scale) in specifications.items():
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"Corrected curvature scale unavailable for {reaction}")
        set_parameter_start(scaled_variables, fitted, scales)
        target = scaled.find_component(target_name)
        target.fix(center * scales[target_name])
        candidates = [center, center - scale, center - 2 * scale, center + scale, center + 2 * scale]
        for candidate in candidates:
            target.set_value(candidate * scales[target_name])
            results = solver.solve(scaled, tee=False)
            if not pyo.check_optimal_termination(results):
                raise RuntimeError(
                    f"Profile failed for {reaction} Delta Cp={candidate}: "
                    f"{results.solver.termination_condition}"
                )
            pyo.TransformationFactory("core.scale_model").propagate_solution(scaled, model)
            evaluation = evaluate(model, fit_blocks, validation)
            objective = pyo.value(model.obj)
            if baseline_objective is None:
                baseline_objective = objective
            external_40_80 = heat_metrics(
                model,
                model.validation_states,
                validation[validation.temperature.isin([40, 80])],
                block_start=0,
            )
            records.append(
                {
                    "reaction": reaction,
                    "dcp_rxn_J_mol_K": candidate,
                    "offset_curvature_scales": (candidate - center) / scale,
                    "objective": objective,
                    "objective_delta": objective - baseline_objective,
                    "vle_n": evaluation["vle"]["n"],
                    "speciation_n": evaluation["speciation"]["n"],
                    **{f"vle_{source}_{metric}": values[metric]
                       for source, values in evaluation["vle"]["sources"].items()
                       for metric in ("n", "pressure_mape_percent", "log_rmse")},
                    **{f"speciation_{source}_{metric}": values[metric]
                       for source, values in evaluation["speciation"]["sources"].items()
                       for metric in ("n", "mae", "rmse", "bias")},
                    **{var.name: pyo.value(var) for var in variables},
                    "vle_pressure_mape_percent": evaluation["vle"]["pressure_mape_percent"],
                    "vle_log_rmse": evaluation["vle"]["log_rmse"],
                    "speciation_mae": evaluation["speciation"]["mae"],
                    "fit_heat_mae": evaluation["fit_heat"]["mae"],
                    "external_40_80_mae": external_40_80["mae"],
                    "external_120_mae": evaluation["external_120"]["mae"],
                    "external_120_rmse": evaluation["external_120"]["rmse"],
                    "external_120_bias": evaluation["external_120"]["bias"],
                }
            )
            print(
                f"{reaction} Delta Cp={candidate:.3f}: "
                f"objective={objective:.4f}, external 120 C MAE="
                f"{evaluation['external_120']['mae']:.4f} kJ/mol"
            )

    result = pd.DataFrame(records).sort_values(["reaction", "dcp_rxn_J_mol_K"])
    result.to_csv(DATA_FILE, index=False)
    return result


def render():
    data = pd.read_csv(DATA_FILE)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    plots = [
        ("external_120_mae", "Kim 2014 120 C MAE (kJ/mol CO2)"),
        ("external_40_80_mae", "Kim 2014 40-80 C MAE (kJ/mol CO2)"),
        ("objective_delta", "Original fit objective change"),
        ("vle_pressure_mape_percent", "Pooled calibration VLE pressure MAPE (%)"),
    ]
    for axis, (column, ylabel) in zip(axes.flat, plots):
        for reaction, group in data.groupby("reaction"):
            axis.plot(
                group.dcp_rxn_J_mol_K,
                group[column],
                marker="o",
                label=reaction,
            )
        axis.set_xlabel("Profiled reaction Delta Cp (J/mol/K)")
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.25)
    axes[0, 0].legend()
    fig.tight_layout()
    fig.savefig(PLOT_FILE, dpi=200)
    plt.close(fig)


def refit_fixed_bicarbonate():
    model, scaled, variables, scaled_variables, fitted, fit_blocks, validation = build_model()
    scales = {
        var.name: iscale.get_scaling_factor(var, default=1) for var in variables
    }
    target_name = "params.reaction_MEA_bicarbonate_formation_combo.dcp_rxn"
    target_value = -284.78887228330444
    set_parameter_start(scaled_variables, fitted, scales)
    target = scaled.find_component(target_name)
    target.fix(target_value * scales[target_name])

    solver = get_solver("ipopt", options=SOLVER_OPTIONS)
    results = solver.solve(scaled, tee=False)
    if not pyo.check_optimal_termination(results):
        raise RuntimeError(
            f"Fixed bicarbonate refit failed: {results.solver.termination_condition}"
        )
    pyo.TransformationFactory("core.scale_model").propagate_solution(scaled, model)

    estimated = [var for var in variables if var.name != target_name]
    estimated_scaled = [scaled.find_component(var.name) for var in estimated]
    estimated_names = {var.name for var in estimated_scaled}
    active_bounds = sum(
        1
        for var, dual in scaled.ipopt_zL_out.items()
        if var.name in estimated_names and dual > 0.2
    ) + sum(
        1
        for var, dual in scaled.ipopt_zU_out.items()
        if var.name in estimated_names and dual < -0.2
    )
    if active_bounds:
        raise RuntimeError(
            f"Fixed bicarbonate refit has {active_bounds} active estimated bounds"
        )
    hessian = np.asarray(get_reduced_hessian(scaled, estimated_scaled), dtype=float)
    eigenvalues, eigenvectors = eigh(hessian)
    if eigenvalues[0] <= 0:
        raise RuntimeError(
            f"Fixed bicarbonate Hessian is not positive definite: {eigenvalues[0]}"
        )
    inverse = eigenvectors @ np.diag(1 / eigenvalues) @ eigenvectors.T
    std_scaled = np.sqrt(np.diag(inverse))
    correlation = inverse / np.outer(std_scaled, std_scaled)
    labels = [var.name for var in estimated]
    pd.DataFrame(correlation, index=labels, columns=labels).to_csv(
        CANDIDATE_CORRELATION_FILE
    )

    curvature = {
        var.name: np.sqrt(inverse[i, i]) / scales[var.name]
        for i, var in enumerate(estimated)
    }
    rows = []
    for var in variables:
        row = fitted[fitted["Object_Name"] == var.name].iloc[0].copy()
        row["Value"] = pyo.value(var)
        row["Curvature_Scale"] = curvature.get(var.name, np.nan)
        row["Relative_Curvature_Scale"] = (
            abs(row["Curvature_Scale"] / row["Value"])
            if np.isfinite(row["Curvature_Scale"])
            else np.nan
        )
        rows.append(row)
    pd.DataFrame(rows, columns=fitted.columns).to_csv(
        CANDIDATE_PARAMETER_FILE, index=False
    )

    evaluation = evaluate(model, fit_blocks, validation)
    summary = pd.DataFrame(
        [
            {
                "candidate": "historical_coordinate_sensitivity",
                "bicarbonate_dcp_J_mol_K": target_value,
                "objective": pyo.value(model.obj),
                "hessian_min_eigenvalue": eigenvalues[0],
                "hessian_max_eigenvalue": eigenvalues[-1],
                "hessian_condition": eigenvalues[-1] / eigenvalues[0],
                "active_estimated_bounds": active_bounds,
                "max_abs_correlation": np.max(
                    np.abs(correlation - np.eye(len(correlation)))
                ),
                "fit_heat_mae": evaluation["fit_heat"]["mae"],
                "external_40_mae": evaluation["external_40"]["mae"],
                "external_80_mae": evaluation["external_80"]["mae"],
                "external_120_mae": evaluation["external_120"]["mae"],
                "vle_pressure_mape_percent": evaluation["vle"]["pressure_mape_percent"],
                "speciation_mae": evaluation["speciation"]["mae"],
            }
        ]
    )
    summary.to_csv(CANDIDATE_SUMMARY_FILE, index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("run", "render", "candidate"))
    args = parser.parse_args()
    {"run": run_profile, "render": render, "candidate": refit_fixed_bicarbonate}[args.stage]()
