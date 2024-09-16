import logging
import pyomo.environ as pyo
import numpy as np

import idaes.logger as idaeslog
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.models_extra.column_models.properties import ModularPropertiesInherentReactionsInitializer
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from Create_Param_Table import create_param_table
from Uncertainty_Analysis import uncertainty_analysis
from Parameter_Setup import get_estimated_params, setup_param_scaling
from eNRTL_property_setup import get_prop_dict
from Load_Datasets import load_datasets
from Plot_Fit import plot_fit

# %% Model Setup

optarg = {
    # 'bound_push' : 1e-22,
    'nlp_scaling_method': 'user-scaling',
    'linear_solver': 'ma57',
    'OF_ma57_automatic_scaling': 'yes',
    'max_iter': 300,
    'tol': 1e-8,
    'constr_viol_tol': 1e-8,
    'halt_on_ampl_error': 'no',
    # 'mu_strategy': 'monotone',
}

param_dic = {'rxn_coeffs': [
    '1',
    '2',
    '3',
    '4',
], 'molecules': [
    'H2O',
    # 'MEA',
    # 'CO2',
], 'cations': [
    'MEAH^+',
], 'anions': [
    'MEACOO^-',
    'HCO3^-',
], 'parameters': [
    'tau_A',
    'tau_B',
    # 'tau_alpha',
], 'interactions': [
    # 'm-ca',
    # 'ca-m',
    # 'm1-m2',
    # 'm2-m1',
    # "ca1-ca2",
    # "ca2-ca1",
],
}
species_dic = {
    'components': ['H2O', 'MEA', 'CO2'],
    'ions': ['MEAH^+', 'MEACOO^-', 'HCO3^-'],
}
system_fit_dic = {
              'temperature': [40.0, 60.0, 80.0, 100.0, 120.0],
              # 'temperature': [40.0],
              'pressure': 101325,
              'amine_weight_percent': .3,
              'loading_constraints': [.1, .6],
              }
column_names = {
    'temperature': 'temperature',
    'loading': 'CO2_loading',
    'amine_concentration': 'MEA_weight_fraction',
    'pressure': 'total_pressure',
    'CO2_pressure': 'CO2_pressure',
    'heat_of_absorption': 'dH_abs'
}
param_table = True


def get_mole_fraction(CO2_loading, amine_concentration):
    MW_MEA = 61.084
    MW_H2O = 18.02

    x_MEA_unloaded = amine_concentration / (MW_MEA / MW_H2O + amine_concentration * (1 - MW_MEA / MW_H2O))
    x_H2O_unloaded = 1 - x_MEA_unloaded

    n_MEA = 100 * x_MEA_unloaded
    n_H2O = 100 * x_H2O_unloaded

    n_CO2 = n_MEA * CO2_loading
    n_tot = n_MEA + n_H2O + n_CO2
    x_CO2, x_MEA, x_H2O = n_CO2 / n_tot, n_MEA / n_tot, n_H2O / n_tot

    mole_frac = {
        'CO2': np.float32(x_CO2),
        'MEA': np.float32(x_MEA),
        'H2O': np.float32(x_H2O),
        'n_T': np.float32(n_tot),
    }
    return mole_frac

if __name__ == "__main__":

    logging.getLogger('pyomo.repn.plugins.nl_writer').setLevel(logging.ERROR)
    init_outlevel = idaeslog.WARNING
    m = pyo.ConcreteModel()
    config = get_prop_dict(species_dic['components'])
    params = m.params = GenericParameterBlock(**config)
    setup_param_scaling(m)

    # %% Dataset Implementation

    obj_expr = 0
    dataset_dir = r"data\data_sets_to_load"

    obj_expr, dfs, param_block_names = load_datasets(m, obj_expr, dataset_dir, species_dic, get_mole_fraction, column_names,
                                                     exclude_list=['Xu', 'Bottinger'])

    # %% Model Initializing and Solving

    iscale.calculate_scaling_factors(m)
    print(f"DOF: {degrees_of_freedom(m)}")
    state_init = ModularPropertiesInherentReactionsInitializer(solver="ipopt",
                                                               solver_options=optarg,
                                                               output_level=init_outlevel)
    for param_block_name in param_block_names:
        param_block = getattr(m, param_block_name)
        state_init.initialize(param_block)

    iscale.calculate_scaling_factors(m)
    df_unfit = get_estimated_params(m, param_dic)
    var_objects = df_unfit['Object'].to_numpy()

    reg = .01
    obj_expr += sum([reg * .5 * ((var - var.value) * iscale.get_scaling_factor(var)) ** 2 for var in var_objects])

    m.obj = pyo.Objective(expr=obj_expr)
    m_scaled = pyo.TransformationFactory('core.scale_model').create_using(m, rename=False)
    m_scaled.ipopt_zL_out = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    m_scaled.ipopt_zU_out = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    var_objects_scaled = [m_scaled.find_component(var.name) for var in var_objects]
    for var in var_objects_scaled:
        var.unfix()
    optarg.pop("nlp_scaling_method", None)  # Scaled model doesn't need user scaling
    solver = get_solver("ipopt", options=optarg)
    solver.solve(m_scaled, tee=False)
    pyo.TransformationFactory('core.scale_model').propagate_solution(m_scaled, m)

    # %% Uncertainty Analysis

    df_fit, W_value = uncertainty_analysis(m_scaled, df_unfit, var_objects, var_objects_scaled)
    if param_table:
        create_param_table(df_fit, title=f"Obj: {pyo.value(m.obj): .4f} - H: {W_value: .4f}")
    obj_value = pyo.value(m.obj)

    # %% Plotting

    plot_fit(df_fit, system_fit_dic, species_dic, get_mole_fraction, obj_value, optarg, config, dataset_dir, column_names)
