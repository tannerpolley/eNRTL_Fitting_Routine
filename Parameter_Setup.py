import idaes.core.util.scaling as iscale
import pandas as pd
import pyomo
import pyomo.environ as pyo
import os


def get_estimated_params(m, fit_param_dic):
    param_dic = {'Description': [], 'Name': [], 'Value': [], 'Uncertainty': [], 'Percent': [], 'Object': [],
                 'Object_Name': []}

    # Pull items that contain the search string
    # Gets all the reactions present in the parameter block as chosen from create_and_scale_params()
    rxns = [item for item in dir(m.params) if "reaction_" in item and "inherent" not in item and "get" not in item]

    obj = getattr(m.params, rxns[0])
    rxn_coeffs_all = [attr for attr in dir(obj) if isinstance(getattr(obj, attr), pyomo.core.base.var.ScalarVar)]

    rxn_coeffs = fit_param_dic['rxn_coeffs']
    molecules = fit_param_dic['molecules']  # List of molecules that to be included in the interactions
    cations = fit_param_dic['cations']  # List of cation to include in the interactions
    anions = fit_param_dic['anions']  # List of anion to include in the interactions
    parameters = fit_param_dic['parameters']  # List of parameters to be used in the fit
    interactions = fit_param_dic['interactions']  # List of interactions accounted for in the fit

    cations_anions = []
    for c in cations:
        for a in anions:
            cations_anions.append(c + ", " + a)

    def latex_mol(mi):
        latex_mi = []
        for char in mi:
            if char.isdigit():
                latex_mi.append('_')
            latex_mi.append(char)
        latex_mi = ''.join(latex_mi)
        return latex_mi

    component_pairs = m.params.Liq.component_pair_set._ordered_values
    component_pairs_chosen = []
    dic = {}

    for m0 in molecules:
        if 'm1-m2' in interactions or 'm2-m1' in interactions:
            beg = '_{'
            end = '}$'
            m1 = m0
            for m2 in molecules:
                if m1 == m2:
                    continue
                else:
                    if 'm1-m2' in interactions:
                        dic[(m1, m2)] = beg + latex_mol(m1) + ',' + latex_mol(m2) + end, 'm1-m2'
                        component_pairs_chosen.append((m1, m2))
                    if 'm2-m1' in interactions:
                        dic[(m2, m1)] = beg + latex_mol(m2) + ',' + latex_mol(m1) + end, 'm2-m1'
                        component_pairs_chosen.append((m2, m1))
        if 'm-ca' in interactions or 'ca-m' in interactions:
            beg = '_{'
            end = ')}$'
            for ca in cations_anions:
                if 'm-ca' in interactions:
                    dic[(m0, ca)] = beg + latex_mol(m0) + ',(' + ca + end, 'm-ca'
                    component_pairs_chosen.append((m0, ca))
                if 'ca-m' in interactions:
                    dic[(ca, m0)] = beg + '(' + ca + '),' + latex_mol(m0) + end, 'ca-m'
                    component_pairs_chosen.append((ca, m0))
    if 'ca1-ca2' in interactions or 'ca2-ca1' in interactions:
        beg = '_{('
        end = ')}$'
        for ca1 in cations_anions:
            for ca2 in cations_anions:
                if ca1 == ca2:
                    continue
                else:
                    if 'ca1-ca2' in interactions:
                        component_pairs_chosen.append((ca1, ca2))
                        dic[(ca1, ca2)] = beg + latex_mol(ca2) + ',' + latex_mol(ca1) + end, 'ca1-ca2'
                    if 'ca2-ca1' in interactions:
                        component_pairs_chosen.append((ca2, ca1))
                        dic[(ca1, ca2)] = beg + latex_mol(ca2) + ',' + latex_mol(ca1) + end, 'ca1-ca2'

    for rxn in rxns:
        rxn_obj = getattr(m.params, rxn)
        rxn_name = rxn_obj.name.split('_')[2][:3]
        for rxn_coeff in rxn_coeffs_all:
            if str(rxn_coeff[-1]) in rxn_coeffs:
                obj = getattr(rxn_obj, rxn_coeff)
                param_dic['Description'].append(f'{rxn_name} rxn coeff')
                param_dic['Object'].append(obj)
                param_dic['Object_Name'].append(obj.name)
                param_dic['Value'].append(pyo.value(obj))
                param_dic['Name'].append('$k_{' + f'{rxn_name},{rxn_coeff[-1]}' + '}$')
                param_dic['Uncertainty'].append(0.0)
                param_dic['Percent'].append(0.0)

    for pair in component_pairs:
        if pair in component_pairs_chosen:
            for parameter in parameters:
                letter = parameter[-1]
                obj = getattr(m.params.Liq, parameter)
                obj = obj[pair]
                param_dic['Description'].append(f"eNRTL {dic[pair][1]}")
                param_dic['Object'].append(obj)
                param_dic['Object_Name'].append(obj.name)
                param_dic['Value'].append(pyo.value(obj))
                param_dic['Name'].append('$' + letter + dic[pair][0])
                param_dic['Uncertainty'].append(0.0)
                param_dic['Percent'].append(0.0)

    folder_path = 'data\Parameters'

    # Check if the folder exists
    if not os.path.exists(folder_path):
        # Create the folder
        os.makedirs(folder_path)

    df_unfit = pd.DataFrame(param_dic)
    df_unfit.to_csv(r'data\Parameters\Parameters_unfit.csv', index=False)

    return df_unfit

def load_fitted_params(m, df):
    # Gets each parameter that was fitted and fixes its value for the run
    for j, row in df.iterrows():
        names = row['Object_Name'].split('.')
        value = row['Value']
        obj = m
        for i in range(len(names) + 1):
            if i == len(names):
                obj.fix(value)
                break
            else:
                name = names[i]
                if name[:3] == 'tau':
                    parameter, species = name[:5], name[5:]
                    obj = getattr(obj, parameter)
                    obj = getattr(obj, '_data')
                    for k, v in obj.items():
                        if v.name == row['Object_Name']:
                            species_key = k
                            break
                    obj = obj[species_key]
                else:
                    obj = getattr(obj, name)


def setup_param_scaling(m):
    params = m.params

    scaling_factor_flow_mol = 1 / 100
    params.set_default_scaling("enth_mol_phase", 3e-4)
    params.set_default_scaling("pressure", 1e-5)
    params.set_default_scaling("temperature", 1)
    params.set_default_scaling("flow_mol", scaling_factor_flow_mol)
    params.set_default_scaling("flow_mol_phase", scaling_factor_flow_mol)
    params.set_default_scaling("flow_mass_phase", scaling_factor_flow_mol / 18e-3)  # MW mixture ~= 24 g/Mol
    params.set_default_scaling("dens_mol_phase", 1 / 18000)
    params.set_default_scaling("visc_d_phase", 700)
    params.set_default_scaling("log_k_eq", 1)

    mole_frac_scaling_factors = {
        "H2O": 2,
        "MEA": 2,
        "CO2": 5
    }
    mole_frac_true_scaling_factors = {
        "MEA": 1e1,
        "H2O": 2,
        "HCO3^-": 5e3,
        "MEAH^+": 1e1,
        "MEACOO^-": 5e2,
        "CO2": 1e3

    }
    inherent_rxn_scaling_factors = {
        "MEA_bicarbonate_formation_combo": 5e2,
        "MEA_carbamate_formation_combo": 5e2,
    }
    for comp, sf_x in mole_frac_scaling_factors.items():
        params.set_default_scaling("mole_frac_comp", sf_x, index=comp)
        params.set_default_scaling("mole_frac_phase_comp", sf_x, index=("Liq", comp))
        params.set_default_scaling("flow_mol_phase_comp", sf_x * scaling_factor_flow_mol, index=("Liq", comp))

    for comp, sf_x in mole_frac_true_scaling_factors.items():
        params.set_default_scaling("mole_frac_phase_comp_true", sf_x, index=("Liq", comp))
        params.set_default_scaling("flow_mol_phase_comp_true", sf_x * scaling_factor_flow_mol, index=("Liq", comp))
    for rxn, sf_xi in inherent_rxn_scaling_factors.items():
        params.set_default_scaling("apparent_inherent_reaction_extent", scaling_factor_flow_mol * sf_xi, index=rxn)

    iscale.set_scaling_factor(m.params.Liq.alpha, 1)
    iscale.set_scaling_factor(m.params.Liq.tau_A, 1)  # Reminder that it's well-scaled
    iscale.set_scaling_factor(m.params.Liq.tau_B, 1 / 300)

    for rxn_name in inherent_rxn_scaling_factors.keys():
        rxn_obj = getattr(m.params, "reaction_" + rxn_name)
        iscale.set_scaling_factor(rxn_obj.k_eq_coeff_1, 1)
        iscale.set_scaling_factor(rxn_obj.k_eq_coeff_2, 1 / 300)
        iscale.set_scaling_factor(rxn_obj.k_eq_coeff_3, 1)
        iscale.set_scaling_factor(rxn_obj.k_eq_coeff_4, 300)
