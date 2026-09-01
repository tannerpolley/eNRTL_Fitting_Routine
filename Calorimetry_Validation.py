import argparse
import os

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
import matplotlib.pyplot as plt
import pandas as pd
import pyomo.environ as pyo
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models_extra.column_models.properties import ModularPropertiesInherentReactionsInitializer

from Fitting_Routine import column_names, get_mole_fraction, species_dic
from Load_Datasets import add_ABS_dataset, heat_of_absorption_expression
from Parameter_Setup import load_fitted_params, setup_param_scaling
from eNRTL_property_setup import get_prop_dict


DATA_FILE = os.path.join("data", "Plots", "Calorimetry_Validation.csv")
PLOT_FILE = os.path.join("data", "Plots", "Calorimetry_Validation.png")


def generate():
    inputs = [
        ("Kim and Svendsen 2007", "kim_2007_dHabs.csv", 0.022),
        ("Kim et al. 2014", "kim_2014_dHabs.csv", None),
    ]
    observations = pd.concat(
        [
            pd.read_csv(os.path.join("data", "data_sets_to_load", filename)).assign(
                source=source,
            )
            for source, filename, _ in inputs
        ],
        ignore_index=True,
    )
    observations["heat_uncertainty_fraction"] = observations.source.map(
        {source: uncertainty for source, _, uncertainty in inputs}
    )
    group_columns = ["source", "temperature", "experiment"]
    groups = list(observations.groupby(group_columns, sort=False))

    model = pyo.ConcreteModel()
    model.params = GenericParameterBlock(**get_prop_dict(species_dic["components"]))
    load_fitted_params(
        model, pd.read_csv(os.path.join("data", "Parameters", "Parameters_fit.csv"))
    )
    setup_param_scaling(model)
    model.states = model.params.build_state_block(
        range(len(observations) + len(groups)), defined_state=True
    )
    add_ABS_dataset(
        model,
        model.states,
        observations,
        column_names,
        species_dic,
        get_mole_fraction,
        0,
    )
    iscale.calculate_scaling_factors(model)
    ModularPropertiesInherentReactionsInitializer(
        solver="ipopt", output_level=idaeslog.WARNING
    ).initialize(model.states)

    records = []
    start = 0
    for (source, _, experiment), group in groups:
        old = model.states[start]
        loading_start = 0.003
        for offset, (_, row) in enumerate(group.iterrows(), start=1):
            current = model.states[start + offset]
            predicted = pyo.value(
                heat_of_absorption_expression(model, old, current)
            ) * 1e-3
            observed = row["dH_abs"]
            records.append(
                {
                    "source": source,
                    "temperature_C": row["temperature"],
                    "experiment": row["experiment"],
                    "role": (
                        "fit"
                        if bool(row["fit"])
                        else "external validation"
                        if source == "Kim et al. 2014"
                        else "temperature holdout"
                    ),
                    "loading_start": loading_start,
                    "loading_end": row["CO2_loading"],
                    "observed_kJ_mol_CO2": observed,
                    "predicted_kJ_mol_CO2": predicted,
                    "residual_kJ_mol_CO2": predicted - observed,
                    "observed_uncertainty_kJ_mol_CO2": (
                        row["heat_uncertainty_fraction"] * observed
                    ),
                }
            )
            old = current
            loading_start = row["CO2_loading"]
        start += len(group) + 1

    result = pd.DataFrame(records)
    assert len(result) == 113
    assert result.groupby(["source", "temperature_C", "experiment"])["loading_end"].apply(
        lambda values: values.is_monotonic_increasing
    ).all()
    result.to_csv(DATA_FILE, index=False)
    for (source, role, temperature), group in result.groupby(
        ["source", "role", "temperature_C"], sort=False
    ):
        print(
            f"{source}, {temperature:g} C, {role}: n={len(group)}, "
            f"MAE={group.residual_kJ_mol_CO2.abs().mean():.4f}, "
            f"RMSE={(group.residual_kJ_mol_CO2.pow(2).mean()) ** 0.5:.4f}, "
            f"bias={group.residual_kJ_mol_CO2.mean():.4f} kJ/mol"
        )


def render():
    data = pd.read_csv(DATA_FILE)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for axis, temperature in zip(axes, (40, 80, 120)):
        subset = data[data.temperature_C == temperature]
        for (source, experiment), group in subset.groupby(["source", "experiment"]):
            color = "C2" if source == "Kim et al. 2014" else f"C{int(experiment) - 1}"
            label = (
                "Kim et al. 2014"
                if source == "Kim et al. 2014"
                else f"Kim-Svendsen 2007 run {int(experiment)}"
            )
            uncertainty = group.observed_uncertainty_kJ_mol_CO2
            if uncertainty.notna().any():
                axis.errorbar(
                    group.loading_end,
                    group.observed_kJ_mol_CO2,
                    yerr=uncertainty,
                    fmt="o",
                    color=color,
                    capsize=2,
                    label=f"{label} observed",
                )
            else:
                axis.scatter(
                    group.loading_end,
                    group.observed_kJ_mol_CO2,
                    marker="o",
                    color=color,
                    label=f"{label} observed",
                )
            axis.scatter(
                group.loading_end,
                group.predicted_kJ_mol_CO2,
                marker="x",
                color=color,
                label=f"{label} model",
            )
        axis.set_title(f"{temperature} °C")
        axis.set_xlabel("CO$_2$ loading (mol/mol MEA)")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("Differential heat of absorption (kJ/mol CO$_2$)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=8)
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    fig.savefig(PLOT_FILE, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("generate", "render"))
    args = parser.parse_args()
    {"generate": generate, "render": render}[args.stage]()
