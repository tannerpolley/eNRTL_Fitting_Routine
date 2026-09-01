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
    observations = pd.read_csv(
        os.path.join("data", "data_sets_to_load", "kim_2007_dHabs.csv")
    )
    groups = list(observations.groupby(["temperature", "experiment"], sort=False))

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
    for (_, experiment), group in groups:
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
                    "temperature_C": row["temperature"],
                    "experiment": row["experiment"],
                    "role": "fit" if bool(row["fit"]) else "holdout",
                    "loading_start": loading_start,
                    "loading_end": row["CO2_loading"],
                    "observed_kJ_mol_CO2": observed,
                    "predicted_kJ_mol_CO2": predicted,
                    "residual_kJ_mol_CO2": predicted - observed,
                    "observed_uncertainty_kJ_mol_CO2": 0.022 * observed,
                }
            )
            old = current
            loading_start = row["CO2_loading"]
        start += len(group) + 1

    result = pd.DataFrame(records)
    assert len(result) == 86
    assert result.groupby(["temperature_C", "experiment"])["loading_end"].apply(
        lambda values: values.is_monotonic_increasing
    ).all()
    result.to_csv(DATA_FILE, index=False)
    for role, group in result.groupby("role"):
        print(
            f"{role}: n={len(group)}, MAE={group.residual_kJ_mol_CO2.abs().mean():.4f}, "
            f"RMSE={(group.residual_kJ_mol_CO2.pow(2).mean()) ** 0.5:.4f}, "
            f"bias={group.residual_kJ_mol_CO2.mean():.4f} kJ/mol"
        )


def render():
    data = pd.read_csv(DATA_FILE)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for axis, temperature in zip(axes, (40, 80, 120)):
        subset = data[data.temperature_C == temperature]
        for experiment, group in subset.groupby("experiment"):
            color = f"C{int(experiment) - 1}"
            axis.errorbar(
                group.loading_end,
                group.observed_kJ_mol_CO2,
                yerr=group.observed_uncertainty_kJ_mol_CO2,
                fmt="o",
                color=color,
                capsize=2,
                label=f"Run {int(experiment)} observed",
            )
            axis.scatter(
                group.loading_end,
                group.predicted_kJ_mol_CO2,
                marker="x",
                color=color,
                label=f"Run {int(experiment)} model",
            )
        role = subset.role.iloc[0]
        axis.set_title(f"{temperature} °C ({role})")
        axis.set_xlabel("CO$_2$ loading (mol/mol MEA)")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("Differential heat of absorption (kJ/mol CO$_2$)")
    axes[-1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(PLOT_FILE, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("generate", "render"))
    args = parser.parse_args()
    {"generate": generate, "render": render}[args.stage]()
