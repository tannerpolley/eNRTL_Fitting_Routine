from scipy.linalg import eigh
import numpy as np
import idaes.core.util.scaling as iscale
import pyomo.environ as pyo
from scipy.sparse.linalg import spsolve
from idaes.core.util.scaling import get_jacobian


def get_reduced_hessian(model, parameters):
    # Assumes a least-squares model has been solved to optimality
    # assert degrees_of_freedom(model) == len(parameters)
    # for val in model.ipopt_zL_out.values():
    #     assert val<1
    # for val in model.ipopt_zU_out.values():
    #     assert val<1
    jac, nlp = get_jacobian(model, scaled=False)
    n_vars = jac.shape[1]
    param_indices = nlp.get_primal_indices(parameters)

    A2 = jac[:, param_indices]

    rest_indices = [i for i in range(n_vars)]

    sorted_indices = param_indices.copy()
    sorted_indices.sort()
    # Go through indices backwards to avoid
    # shifting the list indices when deleting
    # parameter indices
    for idx in sorted_indices[::-1]:
        del rest_indices[idx]

    assert len(rest_indices) == n_vars - len(param_indices)
    for idx in param_indices:
        assert idx not in rest_indices

    A1 = jac[:, rest_indices]
    projection_matrix = spsolve(
        A1.tocsc(),
        A2.tocsc()
    )
    if len(param_indices) > 1:
        # Typically there's extreme fill-in
        # If len(param_indices) == 1, then a (dense) ndarray is returned
        projection_matrix = projection_matrix.todense()
        # import pdb; pdb.set_trace()
    assert nlp.get_obj_factor() == 1

    H = nlp.evaluate_hessian_lag()
    H = H.tocsr()
    H11 = H[np.ix_(rest_indices, rest_indices)]
    H21 = H[np.ix_(param_indices, rest_indices)]
    H12 = H[np.ix_(rest_indices, param_indices)]
    H22 = H[np.ix_(param_indices, param_indices)]

    PT = projection_matrix.transpose()
    H_red = PT @ H11 @ projection_matrix - PT @ H12 - H21 @ projection_matrix + H22
    return H_red


def uncertainty_analysis(m_scaled, df_unfit, estimated_vars, estimated_vars_scaled):
    df_fit = df_unfit.copy()
    active_lbs = {}
    active_ubs = {}
    uncertainty = np.zeros(len(estimated_vars))
    for var, val in m_scaled.ipopt_zL_out.items():
        if val > .2:
            print(var.name, val)
            active_lbs[var.name] = val
            # print('Variables with active lower bounds:')
            # print(var.name, val)
    for val in m_scaled.ipopt_zU_out.values():
        if val < -.2:
            print(var.name, val)
            active_ubs[var.name] = val
            # print('Variables with active upper bounds:')
            # print(var.name, val)
    if len(estimated_vars) > 0 and len(active_lbs) + len(active_ubs) == 0:
        H_red = get_reduced_hessian(m_scaled, estimated_vars_scaled)
        W, V = eigh(H_red)
        inv_red_hess = V @ np.diag(1 / W) @ V.T
        W_value = 1 / W[0]

        if W[0] > 0:
            for i, var in enumerate(estimated_vars):
                uncertainty[i] = pyo.sqrt(inv_red_hess[i][i]) / iscale.get_scaling_factor(var, default=1)
        else:
            print("Warning, Hessian is not positive definite at solution!")
            # print(f"Most negative eigenvalue of reduced Hessian: {W[0]}")
            # print("\n" + "Variables associated:")
            # for i in np.where(abs(V[:, 0]) > 0.3)[0]:
            #     print(str(i) + ": " + estimated_vars_scaled[i].name)
            print("Uncertainties cannot be calculated.")
            for i, var in enumerate(estimated_vars):
                uncertainty[i] = np.NaN
                # print(f"{var.name}: {var.value:.3e}")
    elif len(estimated_vars) > 0:
        print("Warning, active bounds at solution!")
        print("Uncertainties cannot be calculated.")
        for i, var in enumerate(estimated_vars):
            uncertainty[i] = np.NaN
            # print(f"{var.name}: {var.value:.3e}")
        W_value = np.NaN
    else:
        print('Feasibility Problem')
        W_value = np.NaN
        # Feasibility problem
        pass

    for i, var in enumerate(estimated_vars):
        # parameters['Description'].append(var.name)
        # parameters['Name'].append(var.name)
        df_fit.loc[i, 'Value'] = var.value
        df_fit.loc[i, 'Uncertainty'] = uncertainty[i]
        df_fit.loc[i, 'Percent'] = (abs(uncertainty[i] / var.value))
    df_fit.to_csv(r'data\Parameters\Parameters_fit.csv', index=False)

    return df_fit, W_value
