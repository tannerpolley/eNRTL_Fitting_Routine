import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
import pyomo.environ as pyo
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.models_extra.column_models.properties import ModularPropertiesInherentReactionsInitializer
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
import idaes.core.util.scaling as iscale
from Parameter_Setup import setup_param_scaling, load_fitted_params
from eNRTL_property_setup import get_prop_dict
from matplotlib.lines import Line2D


def plot_fit(df, system_fit_dic, species_dic, get_mole_fraction, obj_value, optarg, config, dataset_dir, column_names):
    init_outlevel = idaeslog.WARNING

    def Rochelle_fit(loading, T):
        return np.exp((39.3 - 12155 / T - 19.0 * loading ** 2 + 1105 * loading / T + 12800 * loading ** 2 / T)) / 1e3

    # Assigns a specific marker to a dataset name
    markers = ['o', 's', '^', 'D', 'p', 'h', '*']
    key = {}

    i = 0
    for file in os.listdir(dataset_dir):
        if 'dHabs' not in file:
            name = file.split('_')[0]
            key[name] = markers[i]
            i += 1

    mfc = ['tab:blue', 'tab:orange', 'tab:green', 'tab:cyan', 'tab:purple']
    mfc2 = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']

    Temperature = system_fit_dic['temperature']
    P_sys = system_fit_dic['pressure']
    w_amine = system_fit_dic['amine_weight_percent']
    min_loading_constraint, max_loading_constraint = system_fit_dic['loading_constraints']
    molecules = species_dic['components']
    ions = species_dic['ions']
    molecules_ions = molecules + ions
    molecules_ions.remove('H2O')
    # molecules_ions.remove('CO2')

    temperature = column_names['temperature']
    pressure = column_names['pressure']
    CO2_loading = column_names['loading']
    CO2_pressure = column_names['CO2_pressure']
    amine_concentration = column_names['amine_concentration']
    heat_of_absorption = column_names['heat_of_absorption']


    molecule_key = {}
    for i, molecule in enumerate(molecules_ions):
        molecule_key[molecule] = mfc2[i]

    x_true = {}
    for molecule in molecules_ions:
        x_true[molecule] = []

    lines_model = []
    lines_roch = []
    lines_data = []
    sigma_avg_1 = []
    sigma_avg_2 = []
    fig_VLE, ax_VLE = plt.subplots(figsize=(14, 10))
    fig_Ch_Eq, ax_Ch_Eq = plt.subplots(figsize=(10, 10))

    # Iterates through each temperature value chosen
    for i_t, T in enumerate(Temperature):
        T_K = T + 273.15

        # Gets the min and max loading to set up the loading range for each temperature
        loading_min = []
        loading_max = []
        for file in os.listdir(dataset_dir):
            df_data = pd.read_csv(dataset_dir + '/' + file)
            if 'VLE' in file and (
                    w_amine in df_data[amine_concentration].values and T in df_data[temperature].values):
                df_data = df_data[(df_data[amine_concentration] == w_amine) &
                                  (df_data[temperature] == T) &
                                  (df_data[CO2_loading] > min_loading_constraint) &
                                  (df_data[CO2_loading] < max_loading_constraint)]
                CO2_loading_data = df_data[CO2_loading].to_numpy()
                loading_min.append(min(CO2_loading_data))
                loading_max.append(max(CO2_loading_data))
        loading_constrained = np.linspace(min(loading_min), max(loading_max), 30)
        loading = np.linspace(.01, 1, 30)

        x_true = {}
        for molecule in molecules_ions:
            x_true[molecule] = []

        P_CO2_model = []
        model_data = {
            'P_CO2': [],
            'x_true': x_true
        }
        for alpha in loading:
            m = pyo.ConcreteModel()
            m.params = GenericParameterBlock(**config)
            load_fitted_params(m, df)
            setup_param_scaling(m)
            m.state_block = m.params.build_state_block([0], has_phase_equilibrium=False, defined_state=True)
            blk = m.state_block[0]
            x_dic = get_mole_fraction(alpha, w_amine)
            blk.flow_mol.fix(x_dic['n_T'])
            components = species_dic['components']
            for c in components:
                blk.mole_frac_comp[c].fix(x_dic[c])
            blk.temperature.fix(T_K)
            blk.pressure.fix(P_sys)

            iscale.calculate_scaling_factors(m)
            state_init = ModularPropertiesInherentReactionsInitializer(solver="ipopt",
                                                                       solver_options=optarg,
                                                                       output_level=init_outlevel)
            state_init.initialize(m.state_block)
            m_scaled = pyo.TransformationFactory('core.scale_model').create_using(m, rename=False)
            solver = get_solver("ipopt", options=optarg)
            solver.solve(m_scaled, tee=False)
            pyo.TransformationFactory('core.scale_model').propagate_solution(m_scaled, m)

            model_data['P_CO2'].append(pyo.value(blk.fug_phase_comp["Liq", "CO2"]) / 1e3)
            for molecule in molecules_ions:
                model_data['x_true'][molecule].append(pyo.value(blk.mole_frac_phase_comp_true["Liq", molecule]))

        P_CO2_model_interp = interp1d(loading, model_data['P_CO2'], kind='cubic')
        P_CO2_Roch = Rochelle_fit(loading_constrained, T_K)
        P_CO2_Roch_interp = interp1d(loading_constrained, P_CO2_Roch, kind='cubic')
        data_dir = dataset_dir
        counter = 0

        Sigma_1 = 0
        Sigma_2 = 0
        n = 0
        mec = ['black', 'gray', 'white']
        #%% Plotting data
        for file in os.listdir(data_dir):
            name = file.split('_')[0]
            df_data = pd.read_csv(data_dir + '/' + file)

            if name == 'Bottinger':
                continue

            if name == 'Xu':
                df_data[temperature] = np.round(df_data[temperature].to_numpy(), -1)
                df_data = df_data.sort_values([temperature, CO2_loading])

            #%% VLE data plotting
            if ('VLE' in file
                    and w_amine in df_data[amine_concentration].values
                    and T in df_data[temperature].values):
                df_data = df_data[(df_data[amine_concentration] == w_amine) &
                                  (df_data[temperature] == T) &
                                  (df_data[CO2_loading] > min_loading_constraint) &
                                  (df_data[CO2_loading] < max_loading_constraint)
                                  ]
                loading_data = df_data[CO2_loading].to_numpy()
                CO2_pressure_data = df_data[CO2_pressure].to_numpy()

                for i, alpha in enumerate(loading_data):
                    Sigma_1 += abs((P_CO2_model_interp(alpha) - CO2_pressure_data[i]) / CO2_pressure_data[i])
                    Sigma_2 += abs((P_CO2_Roch_interp(alpha) - CO2_pressure_data[i]) / CO2_pressure_data[i])
                n += len(loading_data)
                data = ax_VLE.semilogy(loading_data, CO2_pressure_data,
                                       label=f"{T} C - {name}", linestyle="none",
                                       marker=key[name], markersize=10, mfc=mfc[i_t], mec=mec[counter])
                lines_data.append(data)

                counter += 1
            #%% Ch Eq data plotting
            if ('ChEq' in file
                    and w_amine in df_data[amine_concentration].values
                    and T in df_data[temperature].values)\
                    and T == 40:
                df_data = df_data[(df_data[amine_concentration] == w_amine) &
                                  (df_data[temperature] == T)
                                  # (df_data[CO2_loading] > min_loading_constraint) &
                                  # (df_data[CO2_loading] < max_loading_constraint)
                                  ]
                loading_data = df_data[CO2_loading].to_numpy()
                molecules_ions = list(df_data.columns)

                molecules_ions.remove(amine_concentration)
                molecules_ions.remove(temperature)
                molecules_ions.remove(CO2_loading)
                # try:
                #     molecules_ions.remove('CO2')
                # except ValueError:
                #     pass
                for molecule in molecules_ions:
                    x_true_i = df_data[molecule].to_numpy()
                    ax_Ch_Eq.plot(loading_data, x_true_i,
                                      label=f"{name}: " + "$x_{" + f"{molecule}" + "}$ data" + f" - T = {T}", linestyle="none",
                                      marker=key[name], markersize=10,
                                      color=molecule_key[molecule])

        model = ax_VLE.semilogy(loading, model_data['P_CO2'], linestyle='dashed', color=mfc[i_t],
                                label=f"{T} C - eNRTL - {Sigma_1 / n:1.2%}")
        sigma_avg_1.append(Sigma_1 / n)
        lines_model.append(model)

        roch = ax_VLE.semilogy(loading_constrained, P_CO2_Roch, linestyle='dotted', color=mfc[i_t],
                               label=f"{T} C - Roch. - {Sigma_2 / n:1.2%}")
        sigma_avg_2.append(Sigma_2 / n)
        lines_roch.append(roch)

        if T == 40:
            for molecule in molecules_ions:
                x_true_i = model_data['x_true'][molecule]
                # loading_ChEq = np.linspace(0, 1, len(x_true_i))
                ax_Ch_Eq.plot(loading, x_true_i,
                                  label="$x_{" + f"{molecule}" + "}$ model", linestyle="dashed",
                                  color=molecule_key[molecule])

    #%% ---- Model Fit

    avg_model = ax_VLE.semilogy([.2], [10], linestyle='-', color='k',
                                label=f"Average MAPE - {np.mean(sigma_avg_1):>2.2%}")
    lines_model.append(avg_model)

    avg_roch = ax_VLE.semilogy([.2], [10], linestyle='-', color='k',
                               label=f"Average MAPE - {np.mean(sigma_avg_2):>2.2%}")
    lines_roch.append(avg_roch)

    handles = [data[0] for data in lines_data]
    labels = [data[0].get_label() for data in lines_data]
    data_legend = ax_VLE.legend(handles, labels, loc='lower right', fontsize='11')
    fig_VLE.gca().add_artist(data_legend)

    #%% ---- Reaction Parameters
    df_rxns = df[df['Description'].str.contains("rxn", case=False, na=False)]
    rxn_descriptions = df_rxns['Description'].to_list()
    rxn_names = [rxn_description.split(" ")[0] for rxn_description in rxn_descriptions]
    rxn_names = list(set(rxn_names))
    df_rxns = [df[df['Description'] == f'{rxn_name} rxn coeff'] for rxn_name in rxn_names]

    add_height = 0
    for df_rxn in df_rxns:
        handles = []
        for i, row in df_rxn.iterrows():
            if abs(row['Value']) < 1:
                name1 = f"{row['Name']} =  $\\bf{{{row['Value']:7.3f}}}$"
            else:
                name1 = f"{row['Name']} =  $\\bf{{{row['Value']:7.1f}}}$"
            handles.append(Line2D([0], [0], label=name1, marker='.', markersize=4,
                                  markerfacecolor='black', linestyle=''))
        parameter_legend = ax_VLE.legend(handles=handles,
                                         loc='upper left', bbox_to_anchor=[1, 1 - add_height], fontsize='10')
        fig_VLE.gca().add_artist(parameter_legend)
        add_height += len(handles) * .033

    #%% --- eNRTL Parameters
    df_eNRTL = df[df['Description'].str.contains("eNRTL", case=False, na=False)]
    handles = []
    for i, row in df_eNRTL.iterrows():
        name = f"{row['Name']} =  $\\bf{{{row['Value']:.3f}}}$"
        handles.append(Line2D([0], [0], label=name, marker='.', markersize=4,
                              markerfacecolor='black', linestyle=''))
    parameter_legend = ax_VLE.legend(handles=handles,
                                     loc='upper left', bbox_to_anchor=[1, 1 - add_height], fontsize='10')
    if len(handles) > 0:
        fig_VLE.gca().add_artist(parameter_legend)

    #%% ---- Rochelle Fit
    handles = [data[0] for data in lines_roch]
    labels = [data[0].get_label() for data in lines_roch]
    data_legend = ax_VLE.legend(handles, labels, loc='upper left', fontsize='11', bbox_to_anchor=[.275, 1])
    fig_VLE.gca().add_artist(data_legend)

    #%% ---- Model Fit
    handles = [model[0] for model in lines_model]
    labels = [model[0].get_label() for model in lines_model]
    ax_VLE.legend(handles, labels, loc='upper left', fontsize='11', )

    #%% Finish Plotting
    ax_VLE.set_xlim(min_loading_constraint, max_loading_constraint)
    ax_VLE.set_ylim(10e-4, 5e3)
    ax_VLE.set_xlabel("CO$_{2}$ Loading, mol CO$_{2}$/mol MEA", fontsize=16)
    ax_VLE.set_ylabel("CO$_{2}$ pressure, kPa", fontsize=16)
    ax_VLE.set_title(
        f"VLE CO2 Solubility at {w_amine:.0%} MEA with {len(df['Value'])} parameters fit - Obj: {obj_value:.2f}",
        fontsize=18)
    ax_VLE.tick_params(labelsize=14)
    fig_VLE.tight_layout()
    fig_VLE.subplots_adjust(right=0.78)

    folder_path = 'data\Plots'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plot_file = "data\Plots\Fitting_Plot.png"
    if os.path.isfile(plot_file):
        os.remove(plot_file)
    fig_VLE.savefig(plot_file)
    print('New plot saved')
    ax_Ch_Eq.set_xlim(0, 1)
    # ax_Ch_Eq.set_ylim(-.05, .125)
    ax_Ch_Eq.set_xlabel("CO$_{2}$ Loading, mol CO$_{2}$/mol MEA", fontsize=16)
    ax_Ch_Eq.set_ylabel("x (mole fraction)", fontsize=16)
    ax_Ch_Eq.set_title(
        f"Speciation at {w_amine:.0%} MEA with {len(df['Value'])} parameters fit - Obj: {obj_value:.2f}",
        fontsize=18)
    ax_Ch_Eq.tick_params(labelsize=14)
    ax_Ch_Eq.legend()
    fig_Ch_Eq.tight_layout()

    plt.show()


if __name__ == '__main__':
    df = pd.read_csv(r'data\Parameters\Parameters_fit.csv')

    from Fitting_Routine import system_fit_dic, species_dic, optarg, column_names, get_mole_fraction

    config = get_prop_dict(["H2O", "MEA", "CO2"])
    dataset_dir = r"data\data_sets_to_load"
    plot_fit(df, system_fit_dic, species_dic, get_mole_fraction, 10.00, optarg, config, dataset_dir, column_names)
