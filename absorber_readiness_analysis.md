# eNRTL fitting routine: formulation and absorber-readiness analysis

Date: 2026-09-01  
Branch: `codex/absorber-readiness-analysis`  
Base: `main` at `b489987` (`fix: make plotting save-only`)

## Executive assessment

**verified:** The reduced two-reaction chemistry, true-species activity equilibrium form, and aqueous-infinite-dilution reference state agree with the Akula MEA-H2O-CO2 formulation. The model equations are suitable for equilibrium/speciation calculations over the fitted domain.

**verified:** The previous four-term equilibrium-constant basis was overparameterized over 313.15-393.15 K. Fixing the unsupported linear-in-temperature term at zero removes two weak coefficient directions while changing the objective by only 0.45% (61.1533 to 61.4272). The new pressure-fit MAPE is 27.51%, compared with about 31.7% from the previous saved result.

**verified:** The eNRTL local-interaction model is incomplete in the default fit. It contains the H2O-MEA molecular parameters, but all molecule-electrolyte terms are held at zero. A bounded test that restored the published fixed interaction values increased the present objective to 1333.36, even after refitting the reaction coefficients. Those terms cannot be switched on independently of the standard-state and calorimetric regression used to obtain them.

**conclusion:** Use the six fitted reaction coefficients as provisional equilibrium initializers. Do not interpret their inverse-Hessian scales as process-design uncertainty, and do not promote the isolated eNRTL interaction fit. A full non-isothermal absorber still requires a validated enthalpy path, transport properties, kinetics, and a vapor package.

## Reproducibility boundary

**verified:** The repository `.venv` contains IDAES 2.7.0.dev0 from the expected fork:

```text
https://github.com/dallan-keylogic/idaes-pse.git
requested revision: eNRTL
commit: 3cdf9b5299b8453fc29274f55d060c1b30a8294b
```

Run with:

```bash
env MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 .venv/bin/python Fitting_Routine.py
```

The complete validated run on this branch took 134.19 s wall time, 131.16 s user CPU time, and 2.58 s system CPU time, with maximum resident memory of 1.46 GB.

## Equation-by-equation review

### Activity basis and reference states

**reference-backed:** Akula et al. define pure-liquid standard states for solvents, infinite-dilution mixed-solvent states for molecular-solute phase equilibrium, and aqueous-infinite-dilution states for chemical equilibrium. The local configuration uses `property_basis="true"`, `ConcentrationForm.activity`, `InfiniteDilutionSingleSolvent`, and H2O as reference component. This is consistent with Akula et al. (2023), p. 3, eqs. 1-3, and Zhang et al. (2011), p. 68, eq. 6.

**verified:** `log_power_law_equil` uses dimensionless true-species activities and `return_log_expression`; therefore the active equilibrium constraints are dimensionally consistent. The unused non-log `return_expression` previously multiplied `exp(logK)` by m3/mol. That dormant inconsistency has been removed so both interfaces now return a dimensionless activity-basis constant.

### Reduced reaction chemistry

**reference-backed:** The implemented reactions are exactly the reduced industrial chemistry in the Akula supporting appendix, p. 18, eqs. D14-D15:

```text
2 MEA + CO2 <-> MEAH+ + MEACOO-
MEA + CO2 + H2O <-> MEAH+ + HCO3-
```

The appendix states that H3O+, OH-, and CO3-- may be neglected at industrially important conditions. Charge, MEA, CO2, and water balances in the code agree with appendix eqs. D18-D21.

**inference:** This reduction is appropriate for the present fitted loading and temperature range. It should not be assumed valid for every stripper or very-low-loading state without a full-chemistry comparison.

### Equilibrium-constant temperature dependence

The retained correlation is:

```text
log(K) = k1 + k2/T + k3*log(T/K)
```

**reference-backed:** Akula derives reaction constants from standard-state Gibbs energies (appendix eq. D22), while the fitted standard-state heat-capacity representation supplies constant heat-capacity terms. That thermodynamic form supports constant, inverse-temperature, and logarithmic-temperature terms. The previous extra `k4*T` term was initialized at zero and was not supported by an additional heat-capacity temperature coefficient in this fit.

**verified:** For the five fitted temperatures, the normalized four-term basis `[1, 1/T, log(T), T]` has condition number `3.208e5`; the three-term basis has condition number `5.813e3`. In the old fit, `k2-k4` correlations were 0.93-0.94 and `k3-k4` correlations were -0.82 to -0.83.

**verified:** Removing `k4` changed the objective from 61.1533 to 61.4272 while reducing the largest relative coefficient scale from 234% to 67%. This is a reduction in unnecessary parameter freedom, not proof that the remaining coefficients are statistically identified.

### eNRTL local interactions

**reference-backed:** The modified symmetric eNRTL formulation contains molecule-molecule, molecule-electrolyte, and electrolyte-electrolyte interactions with `tau = A + B/T`; see Akula et al. (2023), p. 4, Table 1. The local `AkulaTau` rule implements this equation correctly.

**verified:** The default property data provide only H2O-MEA nonzero `tau` terms. Published H2O-(MEAH+, MEACOO-) and H2O-(MEAH+, HCO3-) terms, and their reverse terms, remain commented out. Unspecified interactions default to zero.

**reference-backed:** Akula did not estimate these interaction terms in isolation. The paper simultaneously regressed ion standard-state formation properties, heat-capacity terms, and water-ion-pair interactions against VLE and heat-of-absorption data (pp. 8-11, eq. 34 and Tables 3, 6, and 7). Zhang et al. similarly used binary VLE/excess enthalpy/heat capacity first, then ternary VLE, heat of absorption, heat capacity, and NMR data (2011, pp. 68 and 71-73, Tables 8-10).

**verified:** Adding 12 free eNRTL interaction variables to the old eight reaction variables produced a non-positive-definite reduced Hessian. Restoring the published fixed interactions with only the new six reaction variables produced objective 1333.36. These failures show that the present reaction-only, VLE-plus-speciation objective is not interchangeable with the published joint standard-state regression.

### Enthalpy consistency

**verified:** `enthRxnCullinaneRochelle` is the correct van't Hoff derivative of the implemented log(K):

```text
delta_h_rxn = R*(-k2 + k3*T + k4*T^2)
```

With `k4=0`, this becomes `R*(-k2 + k3*T)`.

**verified:** Liquid enthalpy construction still fails in the installed fork at `idaes/.../eos/enrtl.py:1356`, where a symbolic Pyomo derivative is evaluated in `if dv_dT != 0`. This prevents use of the heat-of-absorption data and blocks a trustworthy non-isothermal column energy balance.

## Objective and data limitations

**verified:** The previous eight-parameter objective contained 240 VLE log-fugacity residuals and 172 speciation residuals. VLE contributed 26.9414, speciation 33.7983, and regularization 0.4137 to total 61.1533.

**verified:** Speciation residuals are multiplied by 10,000, VLE residuals are unweighted log-pressure errors, and no measurement covariance is supplied. The inverse reduced Hessian is therefore a local objective-curvature measure, not a calibrated parameter covariance matrix.

**verified:** `Xu`, `Bottinger`, and `kim` are excluded from fitting. Some excluded VLE data remain visible in the figures. Heat-of-absorption residual construction is present but disabled because of the enthalpy failure.

**inference:** The arbitrary balance between VLE and speciation can shift fitted parameters and makes ordinary standard-error interpretation invalid. Source-level uncertainty weights or a documented sensitivity analysis are needed before formal uncertainty claims.

## Current fit and Hessian

| Reaction | Term | Value | Local curvature scale | Relative |
|---|---|---:|---:|---:|
| bicarbonate | k1 | 176.093 | 9.738 | 5.5% |
| bicarbonate | k2 | -2452.202 | 574.574 | 23.4% |
| bicarbonate | k3 | -27.901 | 1.442 | 5.2% |
| carbamate | k1 | 234.545 | 9.745 | 4.2% |
| carbamate | k2 | -865.310 | 581.988 | 67.3% |
| carbamate | k3 | -37.567 | 1.442 | 3.8% |

**verified:** The new reduced-Hessian eigenvalues span `1.004353e-2` to `5.056428e3`, with condition number `5.034511e5`. The old condition number was `5.205e5`, so the absolute condition improves only 3.3%.

**verified:** The smallest eigenvalue remains at the 0.01 regularization floor. The apparent finite scales for the weakest directions are penalty-influenced. Removing `k4` reduces coefficient ambiguity, but the current data do not create a strongly conditioned statistical estimation problem.

**recommendation:** Propagate uncertainty in predicted `log(K)`, species fractions, or equilibrium pressure at column states, not independent raw-coefficient intervals. Use source/temperature holdouts, profile likelihood, source-cluster bootstrap, or a Bayesian fit only after defensible residual uncertainty models are supplied.

## Fit accuracy and continuation validation

**verified:** The saved prediction grid contains 150 finite states: five temperatures by 30 loadings. The current VLE pressure MAPE values are 39.65%, 34.78%, 17.35%, 26.22%, and 19.53% at 40, 60, 80, 100, and 120 C, respectively; mean MAPE is 27.51%.

**verified:** At 40 C and 30 wt% MEA, the speciation curves retain the same qualitative agreement with Jakobsen NMR data as the old fit. The 20 C speciation data are outside the saved five-temperature grid.

**verified:** Plotting now constructs and initializes one property model per temperature and continues through 30 increasing loadings, rather than constructing and initializing 150 independent models. All 150 solves check optimal termination. Two independently cold-started states agree with the continuation values to maximum relative error `2.674e-9` across CO2 fugacity and six true-species fractions.

**verified:** The complete fit, Hessian, 150-state continuation sweep, table, and two saved figures take 134.19 s. The former complete run took about three minutes, and a fit-only run took about 124 s; most remaining runtime is the simultaneous fit rather than plotting.

## Absorber integration findings

**verified:** The local package does not expose the liquid diffusivity, viscosity, surface tension, and thermal-conductivity methods required by the sibling rate-based MEA column. It also defines no vapor phase.

**verified:** Local ions are named `MEAH^+` and `MEACOO^-`; the sibling absorber eNRTL package uses `MEA_+` and `MEACOO_-`. A deliberate species mapping is required.

**verified:** The active `mea-flowsheet` work uses the same reduced reactions and a custom apparent-enthalpy eNRTL class, but its published molecule-electrolyte terms are also currently commented. Its present tests are construction/property smoke checks, not thermodynamic validation.

**verified:** The MEA-Absorption-Column model currently uses the same six species but concentration-basis equilibrium with a standard-concentration correction. Directly transferring these activity-basis K values would be wrong unless the column uses the same activity convention or performs an explicit standard-state conversion.

**unknown:** Capture, temperature bulge, solvent circulation, pressure drop, and regeneration-duty sensitivity to these fitted constants has not yet been quantified in a solved full column.

## Recommended staged path

1. **Validate the liquid property interface.** Check fugacity, Henry's constant, species fractions, enthalpy, density, viscosity, surface tension, and diffusivities over absorber and stripper states, including units and finite-value checks.
2. **Repair and validate enthalpy.** Correct the symbolic derivative failure, then compare heat capacity and heat of absorption before enabling energy balances or calorimetric regression.
3. **Rebuild the joint thermodynamic regression.** Estimate or fix standard-state ion properties and selected water-ion-pair interactions using VLE plus calorimetry; use speciation as validation or weight it by documented uncertainty.
4. **Align the column interface.** Use one species naming convention, one activity/standard-state convention, validated transport methods, and a separate vapor package.
5. **Bring up the column in stages.** Start with one or two isothermal finite elements, then energy balances, then full discretization and kinetics. Check termination and physical bounds at each stage.
6. **Propagate prediction uncertainty.** Evaluate correlated uncertainty in equilibrium pressure/speciation through capture, rich/lean loading, temperature profile, circulation, and duty only after the thermodynamic fit has defensible uncertainty weights.

## Primary references

- Akula, Lee, Eslick, Bhattacharyya, and Miller, “A Modified Electrolyte Non-Random Two-Liquid Model with Analytical Expression for Excess Enthalpy: Application to the MEA-H2O-CO2 System,” *AIChE Journal* 69(1), 2023, [doi:10.1002/aic.17935](https://doi.org/10.1002/aic.17935). Key locations: pp. 3-4, 8-11; Tables 1, 3, 6, 7; eq. 34.
- Akula et al., “Appendix of eNRTL Model from Akula,” 2023. Key locations: pp. 18-21; eqs. D9-D33.
- Zhang, Que, and Chen, “Thermodynamic Modeling for CO2 Absorption in Aqueous MEA Solution with Electrolyte NRTL Model,” *Fluid Phase Equilibria* 311, 67-75, 2011, [doi:10.1016/j.fluid.2011.08.025](https://doi.org/10.1016/j.fluid.2011.08.025). Key locations: pp. 68, 71-73; Tables 8-10.
- Song and Chen, “Symmetric Electrolyte Nonrandom Two-Liquid Activity Coefficient Model,” *Industrial & Engineering Chemistry Research* 48, 7788-7797, 2009, [doi:10.1021/ie9004578](https://doi.org/10.1021/ie9004578).

## Bottom line

The implemented reduced eNRTL equilibrium formulation is internally sound for its stated activity basis and fitted domain. The six-coefficient fit is simpler and more defensible than the previous eight-coefficient fit, and continuation reduces runtime without changing the selected equilibrium branch. The remaining uncertainty cannot be fixed by adding more free eNRTL parameters: it requires thermodynamically consistent standard-state/calorimetric regression and validated column property methods.
