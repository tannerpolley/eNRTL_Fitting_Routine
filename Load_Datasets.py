import numpy as np
import pyomo.environ as pyo
import os
import pandas as pd
from idaes.models.properties.modular_properties.pure import NIST


def loss(x):
    return 0.5 * x ** 2


def add_VLE_dataset(params, df, column_names, species_dic, get_mole_fraction, obj_expr, has_total_pressure=True):

    for i, row in df.iterrows():
        blk = params[i]
        loading = row[column_names['loading']]
        amine_concentration = row[column_names['amine_concentration']]
        x_dic = get_mole_fraction(loading, amine_concentration)
        blk.flow_mol.fix(x_dic['n_T'])
        components = species_dic['components']
        for c in components:
            blk.mole_frac_comp[c].fix(x_dic[c])
        if has_total_pressure:
            blk.pressure.fix(row[column_names['pressure']] * 1e3)  # Pressure in kPa
        blk.temperature.fix(row[column_names['temperature']] + 273.15)  # Temperature in C
        logfug_CO2 = blk.log_fug_phase_comp["Liq", "CO2"]
        log_P_CO2_data = pyo.log(row[column_names['CO2_pressure']] * 1e3)
        if .1 < loading < .6:
            obj_expr += loss(logfug_CO2 - log_P_CO2_data)  # Pressure in kPa))

    return obj_expr


def add_ABS_dataset(m, params, df, column_names, species_dic, get_mole_fraction, obj_expr):
    unique_temperatures = df[column_names['temperature']].unique()
    dfs = [df[df[column_names['temperature']] == t] for t in unique_temperatures]

    idx_start = 0
    for count, df in enumerate(dfs):

        T = df[column_names['temperature']].unique()[0] + 273.15

        # Start out with a completely unloaded mixture
        blk = params[idx_start]
        amine_concentration = df.iloc[0][column_names['amine_concentration']]
        x_dic = get_mole_fraction(.003, amine_concentration)
        blk.flow_mol.fix(x_dic['n_T'])
        components = species_dic['components']
        for c in components:
            blk.mole_frac_comp[c].fix(x_dic[c])
        # Not a lot of information about pressure
        if 'pressure' in list(df.columns):
            blk.pressure.fix(df.iloc[0][column_names['pressure']])
        else:
            blk.pressure.fix(101325)
        blk.temperature.fix(T)

        for i, row in df.iterrows():
            blk_old = blk
            blk = params[count + i + 1]
            loading = row[column_names['loading']]
            amine_concentration = row[column_names['amine_concentration']]
            x_dic = get_mole_fraction(loading, amine_concentration)
            blk.flow_mol.fix(x_dic['n_T'])
            for c in components:
                blk.mole_frac_comp[c].fix(x_dic[c])
            # Not a lot of information about pressure
            if 'pressure' in list(df.columns):
                blk.pressure.fix(row[column_names['pressure']])
            else:
                blk.pressure.fix(101325)
            blk.temperature.fix(T)
            CO2_obj = m.params.CO2
            Hl_f = blk.energy_internal_mol_phase["Liq"]
            F_f = blk.flow_mol
            Hl_i = blk_old.energy_internal_mol_phase["Liq"]
            F_i = blk_old.flow_mol
            H_ig = NIST.enth_mol_ig_comp.return_expression(blk, CO2_obj, blk.temperature)

            Ff_CO2 = blk.flow_mol_comp["CO2"]
            Fi_CO2 = blk_old.flow_mol_comp["CO2"]
            dH_abs_expr = -(Hl_f * F_f - Hl_i * F_i - H_ig * (Ff_CO2 - Fi_CO2)) / (Ff_CO2 - Fi_CO2)

            # if (row['CO2_loading'] <= 0.4) and (row[column_names['heat_of_absorption']] <= 130):  # threshold bc of missing vapor phase enthalpy
            #     residual_scale = 50
            #     residual = (dH_abs_expr*1e-3 - row[column_names['heat_of_absorption']])/residual_scale
            #     obj_expr += loss(residual)

        idx_start += len(df[column_names['loading']]) + 1
    return obj_expr


def add_ChEq_dataset(params, df, column_names, species_dic, get_mole_fraction, obj_expr, skip_CO2_speciation):
    for i, row in df.iterrows():
        blk = params[i]
        T = row[column_names['temperature']]
        loading = row[column_names['loading']]
        amine_concentration = row[column_names['amine_concentration']]
        x_dic = get_mole_fraction(loading, amine_concentration)
        blk.flow_mol.fix(x_dic['n_T'])
        blk.temperature.fix(T + 273.15)  # Temperature in C
        components = species_dic['components']
        for c in components:
            blk.mole_frac_comp[c].fix(x_dic[c])

        molecules_ions = list(df.columns)
        molecules_ions.remove(column_names['amine_concentration'])
        molecules_ions.remove(column_names['temperature'])
        molecules_ions.remove(column_names['loading'])
        if skip_CO2_speciation:
            try:
                molecules_ions.remove('CO2')
            except ValueError:
                pass

        for molecule in molecules_ions:
            x_true_model = blk.mole_frac_phase_comp_true["Liq", molecule]
            x_true_data = row[molecule]
            obj_expr += loss(x_true_model - x_true_data)*10000
    return obj_expr


def load_datasets(m, obj_expr, dataset_dir, species_dic, get_mole_fraction, column_names, skip_CO2_speciation=True, exclude_list=None):
    if exclude_list is None:
        exclude_list = []
    dfs = []
    param_block_names = []
    for name in os.listdir(dataset_dir):
        filename = dataset_dir + os.sep + name
        name, year, dataset_type = name.split('_')
        if name in exclude_list:
            continue
        dataset_type = dataset_type.split('.')[0]
        df = pd.read_csv(filename, index_col=None)
        dfs.append(df)
        param_block_name = name + '_' + dataset_type
        param_block_names.append(param_block_name)

        if dataset_type == 'dHabs':
            setattr(m, param_block_name, m.params.build_state_block(range(len(df) + len(df[column_names['temperature']].unique())),
                                                                    defined_state=True))
        else:
            setattr(m, param_block_name, m.params.build_state_block(range(len(df)), defined_state=True))

        param_block = getattr(m, param_block_name)
        if dataset_type == 'VLE':
            if 'total_pressure' in df.columns:
                obj_expr = add_VLE_dataset(param_block, df, column_names, species_dic, get_mole_fraction, obj_expr)
            else:
                obj_expr = add_VLE_dataset(param_block, df, column_names, species_dic, get_mole_fraction, obj_expr, has_total_pressure=False)
        elif dataset_type == 'dHabs':
            obj_expr = add_ABS_dataset(m, param_block, df, column_names, species_dic, get_mole_fraction, obj_expr)
        elif dataset_type == 'ChEq':
            obj_expr = add_ChEq_dataset(param_block, df, column_names, species_dic, get_mole_fraction, obj_expr, skip_CO2_speciation)

    return obj_expr, dfs, param_block_names
