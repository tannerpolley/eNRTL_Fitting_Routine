"""Fixed-true-composition Debye–Hückel invariant; no equilibrium re-solve."""
import unittest

import numpy as np
import pyomo.environ as pyo
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.properties.modular_properties.eos.enrtl import ClosestApproach, ENRTL
from pyomo.core.expr.calculus.derivatives import Modes, differentiate

from eNRTL_property_setup import ENRTLDirectTemperatureDerivative, get_prop_dict


def derivative_sweep():
    model = pyo.ConcreteModel()
    model.params = GenericParameterBlock(**get_prop_dict(["H2O", "MEA", "CO2"]))
    model.states = model.params.build_state_block([0], defined_state=True)
    state = model.states[0]
    composition = {"H2O": .85, "MEA": .09, "CO2": .01,
                   "MEAH^+": .025, "MEACOO^-": .015, "HCO3^-": .01}
    for species, fraction in composition.items():
        state.mole_frac_phase_comp_true["Liq", species].fix(fraction)
    A = state.Liq_A_DH
    v = state.Liq_vol_mol_solvent
    eps = state.Liq_relative_permittivity_solvent
    old = -A / 2 * (
        differentiate(v, state.temperature, mode=Modes.reverse_symbolic) / v
        + 3 * state.Liq_d_relative_permittivity_solvent_dT / eps
    )
    analytic = ENRTLDirectTemperatureDerivative._dA_DH_dT(state, "Liq")
    records = []
    for temperature in (313.15, 353.15, 393.15):
        state.temperature.set_value(temperature)
        corrected, historical = pyo.value(analytic), pyo.value(old)
        for step in (1., .1, .01, .001, .0001):
            state.temperature.set_value(temperature + step)
            plus = pyo.value(A)
            state.temperature.set_value(temperature - step)
            minus = pyo.value(A)
            state.temperature.set_value(temperature)
            finite_difference = (plus - minus) / (2 * step)
            records.append(dict(temperature_K=temperature, step_K=step,
                                central_difference_per_K=finite_difference,
                                corrected_per_K=corrected, old_per_K=historical))
    return records


def enthalpy_debye_huckel_sweep():
    model = pyo.ConcreteModel()
    model.params = GenericParameterBlock(**get_prop_dict(["H2O", "MEA", "CO2"]))
    model.states = model.params.build_state_block([0], defined_state=True)
    state = model.states[0]
    for species, fraction in {
        "H2O": 0.85,
        "MEA": 0.09,
        "CO2": 0.01,
        "MEAH^+": 0.025,
        "MEACOO^-": 0.015,
        "HCO3^-": 0.01,
    }.items():
        state.mole_frac_phase_comp_true["Liq", species].fix(fraction)
    R = pyo.value(ENRTL.gas_constant(state))
    ionic_strength = state.Liq_ionic_strength
    factor = pyo.value(
        4
        * ionic_strength
        / ClosestApproach
        * pyo.log(1 + ClosestApproach * ionic_strength**0.5)
    )
    A = state.Liq_A_DH
    records = []
    for temperature in (313.15, 353.15, 393.15):
        state.temperature.set_value(temperature)
        routed = pyo.value(
            R
            * temperature**2
            * factor
            * ENRTLDirectTemperatureDerivative._dA_DH_dT(state, "Liq")
        )
        total_enthalpy = pyo.value(state.enth_mol_phase["Liq"])
        for step in (1.0, 0.01):
            state.temperature.set_value(temperature + step)
            plus = pyo.value(A)
            state.temperature.set_value(temperature - step)
            minus = pyo.value(A)
            state.temperature.set_value(temperature)
            finite_difference = R * temperature**2 * factor * (plus - minus) / (2 * step)
            records.append(
                {
                    "temperature_K": temperature,
                    "step_K": step,
                    "finite_difference": finite_difference,
                    "routed": routed,
                    "total_enthalpy": total_enthalpy,
                }
            )
    return records


class DebyeHuckelInvariant(unittest.TestCase):
    def test_step_sweep(self):
        for row in derivative_sweep():
            fd = row["central_difference_per_K"]
            # Central differences have O(step**2) truncation error; at 1 K
            # 1e-4 relative agreement suffices, with 1e-7 on the fine plateau.
            tolerance = 1e-4 if row["step_K"] == 1 else 1e-6
            self.assertLess(abs(row["corrected_per_K"] / fd - 1), tolerance)
            self.assertGreater(abs(row["old_per_K"] / fd - 1), .1)

    def test_production_excess_enthalpy_route(self):
        for row in enthalpy_debye_huckel_sweep():
            tolerance = 1e-4 if row["step_K"] == 1 else 1e-6
            self.assertLess(abs(row["routed"] / row["finite_difference"] - 1), tolerance)
            self.assertTrue(np.isfinite(row["total_enthalpy"]))


if __name__ == "__main__":
    unittest.main()
