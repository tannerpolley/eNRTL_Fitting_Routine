"""Run a small MEA-column check using fitted calorimetry parameters.

The column implementation lives in the sibling ``mea-flowsheet`` checkout.
This adapter changes only an in-memory property-package configuration and
writes validation results in this repository.
"""

import argparse
import contextlib
import io
import math
import os
import sys
import time
from pathlib import Path

import pandas as pd


R = 8.31446261815324
T_REF = 353.15
DEFAULT_PARAMETER_FILE = Path("data/Parameters/Parameters_fixed_bicarbonate.csv")
DEFAULT_SUMMARY_FILE = Path("data/Plots/Absorber_Calorimetry_Summary.csv")


def physical_to_legacy_coefficients(parameters):
    """Convert (ln K_ref, Delta H_ref, Delta Cp) to the column basis."""
    result = {}
    for reaction, values in parameters.items():
        log_k_ref, dh_ref, dcp = values
        k3 = dcp / R
        k2 = (-dh_ref + dcp * T_REF) / R
        k1 = log_k_ref - k2 / T_REF - k3 * math.log(T_REF)
        result[reaction] = (k1, k2, k3, 0.0)
    return result


def load_parameters(path):
    table = pd.read_csv(path)
    values = {}
    for reaction in ("bicarbonate", "carbamate"):
        rows = table[table["Object_Name"].str.contains(reaction)]
        if rows.empty:
            raise ValueError(f"No {reaction} reaction in {path}")
        values[f"MEA_{reaction}_formation_combo"] = (
            float(rows[rows["Object_Name"].str.endswith("log_k_ref")]["Value"].iloc[0]),
            float(rows[rows["Object_Name"].str.endswith("dh_rxn_ref")]["Value"].iloc[0]),
            float(rows[rows["Object_Name"].str.endswith("dcp_rxn")]["Value"].iloc[0]),
        )
    return values


def import_column(mea_root, idaes_root, legacy_coefficients):
    """Import the sibling column and inject coefficients before construction."""
    sys.path[:0] = [str(idaes_root), str(mea_root / "flowsheets"), str(mea_root)]
    os.chdir(mea_root)

    import mea_properties

    original_generic_package = mea_properties.GenericParameterBlock

    def configured_generic_package(**config):
        reactions = config.get("inherent_reactions", {})
        for reaction, coefficients in legacy_coefficients.items():
            if reaction not in reactions:
                continue
            entries = reactions[reaction]["parameter_data"]["k_eq_coeff"]
            for index, coefficient in enumerate(coefficients, start=1):
                key = str(index)
                entries[key] = (coefficient, entries[key][1])
        return original_generic_package(**config)

    mea_properties.GenericParameterBlock = configured_generic_package

    import MEA_solvent_absorber_model_NGCC_RT_F_CO2loading as column

    return column


def configure_model(column, temperature, nfe):
    import pyomo.environ as pyo
    import idaes.core.util.scaling as iscale
    from idaes.core.util.model_statistics import unused_variables_set
    from eNRTL.MEA_eNRTL import initialize_inherent_reactions

    model = column.build_column_model(column.get_uniform_grid(nfe), False, True)
    absorber = model.fs.absorber
    absorber.diameter_column.fix(12)
    absorber.length_column.fix(20)
    absorber.vapor_inlet.flow_mol.fix(12000)
    absorber.vapor_inlet.temperature.fix(temperature)
    absorber.vapor_inlet.pressure.fix(105000)
    for component, fraction in {
        "CO2": 0.042,
        "H2O": 0.058,
        "N2": 0.77,
        "O2": 0.13,
    }.items():
        absorber.vapor_inlet.mole_frac_comp[0, component].fix(fraction)

    co2_loading = 0.20
    water_loading = 7.878958
    total = co2_loading + water_loading + 1
    absorber.liquid_inlet.flow_mol.fix(25000)
    absorber.liquid_inlet.temperature.fix(temperature)
    absorber.liquid_inlet.pressure.fix(258700)
    for component, fraction in {
        "MEA": 1 / total,
        "CO2": co2_loading / total,
        "H2O": water_loading / total,
    }.items():
        absorber.liquid_inlet.mole_frac_comp[0, component].fix(fraction)

    column.define_column_design_parameters(model)
    column.scale_mea_liquid_params(
        model.fs.liquid_properties, ions=True, scaling_factor_flow_mol=3e-4
    )
    column.scale_mea_vapor_params(
        model.fs.vapor_properties, scaling_factor_flow_mol=3e-4
    )

    saved_bounds = [
        (var, var.lb, var.ub, var.domain)
        for var in model.component_data_objects(pyo.Var, descend_into=True)
    ]
    absorber.vapor_phase.estimate_states(always_estimate=True)
    absorber.liquid_phase.estimate_states(always_estimate=True)
    initialize_inherent_reactions(absorber.liquid_phase.properties)
    column.strip_statevar_bounds_for_initialization(model)
    iscale.calculate_scaling_factors(model)
    transform = pyo.TransformationFactory("contrib.strip_var_bounds")
    transform.apply_to(model, reversible=True)
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
            io.StringIO()
        ):
            absorber.initialize(
                outlvl=0,
                optarg={"linear_solver": "ma57", "max_iter": 1200, "tol": 1e-8},
            )
    finally:
        transform.revert(model)
        for var, lower, upper, domain in saved_bounds:
            var.setlb(lower)
            var.setub(upper)
            var.domain = domain
    for var in unused_variables_set(model):
        var.fix()
    return model


def solve_model(model, solver, label):
    import pyomo.environ as pyo
    from idaes.core.util.model_statistics import degrees_of_freedom

    start = time.perf_counter()
    result = solver.solve(model, tee=False)
    if not pyo.check_optimal_termination(result):
        raise RuntimeError(f"{label} solve failed: {result.solver.termination_condition}")

    absorber = model.fs.absorber
    vapor_in = absorber.vapor_phase.properties[0, 0]
    vapor_out = absorber.vapor_phase.properties[0, 1]
    liquid_in = absorber.liquid_phase.properties[0, 1]
    liquid_out = absorber.liquid_phase.properties[0, 0]
    max_residual = 0.0
    for constraint in model.component_data_objects(
        pyo.Constraint, active=True, descend_into=True
    ):
        if constraint.has_lb():
            max_residual = max(
                max_residual, abs(pyo.value(constraint.body - constraint.lower))
            )
        if constraint.has_ub():
            max_residual = max(
                max_residual, abs(pyo.value(constraint.body - constraint.upper))
            )
    mole_fractions = [
        pyo.value(var)
        for var in model.component_data_objects(pyo.Var, descend_into=True)
        if "mole_frac" in var.name and var.value is not None
    ]

    return {
        "temperature_C": pyo.value(absorber.vapor_inlet.temperature[0]) - 273.15,
        "label": label,
        "termination": str(result.solver.termination_condition),
        "solve_s": time.perf_counter() - start,
        "degrees_of_freedom": degrees_of_freedom(model),
        "capture_percent": pyo.value(absorber.co2_capture[0]),
        "co2_vapor_out_mol_s": pyo.value(
            vapor_out.flow_mol_phase_comp["Vap", "CO2"]
        ),
        "liquid_outlet_temperature_C": pyo.value(liquid_out.temperature) - 273.15,
        "vapor_outlet_temperature_C": pyo.value(vapor_out.temperature) - 273.15,
        "liquid_inlet_temperature_C": pyo.value(liquid_in.temperature) - 273.15,
        "vapor_inlet_temperature_C": pyo.value(vapor_in.temperature) - 273.15,
        "max_constraint_residual": max_residual,
        "minimum_mole_fraction": min(mole_fractions),
        "co2_vapor_in_mol_s": pyo.value(
            vapor_in.flow_mol_phase_comp["Vap", "CO2"]
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mea-flowsheet-root", type=Path, default=Path("../mea-flowsheet")
    )
    parser.add_argument("--idaes-pse-root", type=Path, default=Path("../idaes-pse"))
    parser.add_argument("--parameter-file", type=Path, default=DEFAULT_PARAMETER_FILE)
    parser.add_argument("--summary-file", type=Path, default=DEFAULT_SUMMARY_FILE)
    parser.add_argument(
        "--temperatures-K", nargs="+", type=float, default=[313.15, 353.15]
    )
    parser.add_argument("--nfe", type=int, default=40)
    args = parser.parse_args()

    fitting_root = Path(__file__).resolve().parent
    parameter_file = (fitting_root / args.parameter_file).resolve()
    summary_file = (fitting_root / args.summary_file).resolve()
    mea_root = (fitting_root / args.mea_flowsheet_root).resolve()
    idaes_root = (fitting_root / args.idaes_pse_root).resolve()

    physical = load_parameters(parameter_file)
    legacy = physical_to_legacy_coefficients(physical)
    column = import_column(mea_root, idaes_root, legacy)
    from idaes.core.solvers import get_solver

    coefficient_file = fitting_root / (
        "data/Parameters/Parameters_fixed_bicarbonate_legacy_coefficients.csv"
    )
    rows = []
    for reaction, coefficients in legacy.items():
        for index, coefficient in enumerate(coefficients, start=1):
            rows.append(
                {"reaction": reaction, "coefficient": index, "value": coefficient}
            )
    coefficient_file.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(coefficient_file, index=False)

    model = configure_model(column, args.temperatures_K[0], args.nfe)
    solver = get_solver(
        "ipopt",
        options={
            "nlp_scaling_method": "user-scaling",
            "linear_solver": "ma57",
            "OF_ma57_automatic_scaling": "yes",
            "max_iter": 1200,
            "tol": 1e-8,
            "constr_viol_tol": 1e-8,
            "halt_on_ampl_error": "no",
        },
    )
    records = [solve_model(model, solver, "cold_start")]
    for temperature in args.temperatures_K[1:]:
        model.fs.absorber.vapor_inlet.temperature[0].set_value(temperature)
        model.fs.absorber.liquid_inlet.temperature[0].set_value(temperature)
        records.append(solve_model(model, solver, "warm_start"))

    summary_file.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(summary_file, index=False)
    print(pd.DataFrame(records).to_string(index=False))
    print(f"legacy coefficients: {coefficient_file}")
    print(f"absorber summary: {summary_file}")


if __name__ == "__main__":
    main()
