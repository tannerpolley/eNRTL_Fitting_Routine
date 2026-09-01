# eNRTL fitting routine: absorber-readiness analysis

Date: 2026-09-01  
Branch: `codex/absorber-readiness-analysis`  
Base: `main` at `b489987` (`fix: make plotting save-only`)

## Executive assessment

**Verified:** The current routine is usable as an equilibrium/speciation regression and produces the requested CO2-pressure and speciation figures. The default run has eight fitted reaction-equilibrium correlation coefficients, 412 active residual terms, zero degrees of freedom before optimization, and a final objective of approximately 61.1533.

**Verified:** The fitted reaction coefficients are not kinetic rate constants. They define temperature-dependent equilibrium constants for the bicarbonate and carbamate reaction combinations. A rate-based absorber still needs independently validated reaction kinetics, liquid/vapor mass-transfer correlations, packing/hydrodynamics, and a vapor property package.

**Conclusion:** These parameters are reasonable initial equilibrium parameters for a staged absorber bring-up, but they should not yet be treated as validated parameters for a non-isothermal, rate-based absorber/stripper or for process design. The strongest current blockers are the liquid enthalpy construction failure and missing local transport-property interfaces.

## Reproducibility boundary

**Verified:** The fitting environment is the repository `.venv`, with IDAES installed from the expected fork and branch:

```text
https://github.com/dallan-keylogic/idaes-pse.git
requested revision: eNRTL
commit: 3cdf9b5299b8453fc29274f55d060c1b30a8294b
version: 2.7.0.dev0
```

The default command is:

```bash
env MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 .venv/bin/python Fitting_Routine.py
```

The local analysis branch started clean from `main`; `main` and the sibling `idaes-pse` checkout were not modified.

## What is currently fitted

The default `Fitting_Routine.py` enables only the four reaction-correlation terms for each of two reaction combinations:

```text
log(K) = k1 + k2/T + k3*log(T/K) + k4*T
```

The eNRTL interaction parameter lists are empty by default. Therefore the default results fit reaction-equilibrium coefficients against VLE and speciation data while holding the eNRTL interaction parameters fixed at their configured values.

**Verified:** The active objective contains:

| Contribution | Terms | Objective contribution | Share |
|---|---:|---:|---:|
| VLE log-fugacity residuals | 240 | 26.9414 | 44.1% |
| Speciation residuals | 172 | 33.7983 | 55.3% |
| Quadratic regularization | 8 parameters | 0.4137 | 0.7% |
| Total | 412 | 61.1533 | 100% |

The speciation residuals are multiplied by 10,000, while the VLE residuals are log-pressure residuals. No measurement standard deviations are used. Consequently, the objective is a useful optimization score, but it is not yet a statistically interpretable chi-square or likelihood.

**Verified:** The current fit excludes `Xu`, `Bottinger`, and `kim` from the objective. `kim` is excluded because constructing its heat-of-absorption state triggers the installed eNRTL/Pyomo symbolic-Boolean failure. The plots can still display some excluded data, so plotted coverage and fitted coverage are not identical.

## Fit accuracy and domain coverage

**Verified:** The saved model grid contains 150 finite states (five temperatures by 30 loadings). The six stored true-species mole fractions sum to one within `8.6e-10`, and no negative species fractions were observed.

**Verified:** On an independent recomputation using the saved grid and the 30 wt% VLE observations within the plot range, the model has approximately 31.7% loading-weighted pressure MAPE and RMS log error 0.351. Per-temperature pressure MAPE is approximately 40.2%, 36.4%, 17.3%, 28.0%, and 21.1% at 40, 60, 80, 100, and 120 °C, respectively. This is a useful screening-level fit, not a high-accuracy validation result.

**Verified:** At 40 °C, 30 wt% MEA, interpolation against the Jakobsen speciation data gives mean absolute errors of approximately 0.004–0.008 mole fraction for the fitted species, with maximum absolute errors of approximately 0.012–0.017. The 20 °C speciation data are not represented in the saved five-temperature prediction grid, so low-temperature behavior remains less checked.

**Inference:** The fit is strongest as an interpolation model over the data-supported equilibrium region. Extrapolation to absorber/stripper temperature and loading profiles should be checked explicitly, especially at low loading and at temperatures where only a small number of sources contribute.

## Parameter uncertainty and identifiability

The current uncertainty calculation is a local, linearized reduced-Hessian calculation. It is not a bootstrap, profile likelihood, Bayesian posterior, or global uncertainty analysis.

**Verified:** The default fitted values and reported one-sigma-like local uncertainties are:

| Reaction | Term | Value | Local uncertainty | Relative |
|---|---|---:|---:|---:|
| bicarbonate | k1 | 176.675 | 9.788 | 5.5% |
| bicarbonate | k2 | -1534.853 | 1568.784 | 102.2% |
| bicarbonate | k3 | -29.186 | 2.490 | 8.5% |
| bicarbonate | k4 | 0.012316 | 0.019078 | 154.9% |
| carbamate | k1 | 234.214 | 9.789 | 4.2% |
| carbamate | k2 | -1515.788 | 1705.344 | 112.5% |
| carbamate | k3 | -36.665 | 2.609 | 7.1% |
| carbamate | k4 | -0.008835 | 0.020674 | 234.0% |

**Verified:** The reduced-Hessian eigenvalues span approximately `1.0e-2` to `5.2e3`, for an absolute condition number of approximately `5.2e5`. The smallest two eigenvalues are at the regularization floor of about `0.01`.

**Inference:** The `k2` and `k4` terms are weakly identified and should not be interpreted individually. The largest parameter correlations include `k2`–`k4` of about `+0.93` to `+0.94`, and `k3`–`k4` of about `-0.82` to `-0.83`. The temperature basis is allowing multiple coefficient combinations to produce similar `K(T)` values over the measured range.

For absorber use, uncertainty should therefore be propagated through predicted `log(K)` or equilibrium pressure at the actual column states, rather than through independent uncertainty intervals on the eight raw coefficients. Profile likelihoods, dataset/temperature holdouts, multistart fits, bootstrap resampling by source, or a Bayesian posterior would be more defensible before reporting process-level uncertainty.

**Verified:** An isolated eNRTL-on smoke fit completed after adding 12 eNRTL variables to the eight reaction variables, but its reduced Hessian was not positive definite and the uncertainty routine returned `NaN` uncertainties. The eNRTL-on parameter values from that smoke fit should not be promoted to validated parameters.

## Efficiency findings

**Verified:** The fit and uncertainty stage took approximately 124 seconds in a fresh no-plot diagnostic run. A complete default run took approximately three minutes.

**Verified:** `Plot_Fit.py` performs 5 temperature sweeps × 30 loading points = 150 fresh property-model constructions and IPOPT solves. Each point creates a new `GenericParameterBlock`, state block, initializer, scaled model, and solver call. The previous point's solution is not used as a continuation start, and solver termination is not checked before values are plotted.

**Recommendation:** The highest-value efficiency improvement is to retain one model per temperature and warm-start successive loadings, or at minimum reuse the previous point's state values as initialization. A second useful improvement is to separate fitting from plotting so the saved fitted parameters and comparison CSV can be replotted without repeating 150 solves. Both changes should preserve the current numerical results before being adopted.

## Full absorber readiness

### Hard blockers

**Verified:** Building the local liquid state and requesting `enth_mol_phase` fails in the installed fork at `idaes/.../eos/enrtl.py:1356`, where a symbolic Pyomo derivative is evaluated in a Python Boolean test:

```python
if dv_dT != 0:
```

The resulting `PyomoException` prevents enthalpy construction. This is also why the heat-of-absorption dataset is currently excluded. A non-isothermal absorber, and especially a stripper energy balance, cannot be trusted until this path is corrected and validated.

**Verified:** The local fitting configuration comments out phase transport hooks and component diffusivity hooks. A state-block probe could not retrieve `diffus_phase_comp["Liq", "CO2"]`, `diffus_phase_comp["Liq", "MEA"]`, or the ionic diffusivities, and `visc_d_phase`, `surf_tens_phase`, and `therm_cond_phase` were not provided by the package configuration.

**Verified:** The local package defines only a liquid phase. The sibling `MEAColumn` implementation uses separate vapor and liquid property packages and directly requests liquid surface tension, viscosity, and diffusivity for its rate-based mass-transfer equations.

**Verified:** The local ion names are `MEAH^+` and `MEACOO^-`, while the sibling absorber implementation expects names such as `MEA_+` and `MEACOO_-`. This is an integration contract mismatch, not just a plotting-label difference.

### Important modeling limits

**Inference:** The sibling `ENRTLApparentEnthalpy` class is a possible starting point for the enthalpy repair, but it changes the enthalpy route to an apparent-species sum. It must be checked against heat capacity, heat of absorption, energy balances, and thermodynamic consistency before being used in a process model; copying it without those checks would only hide the construction failure.

**Verified:** The current routine does not fit reaction kinetics, mass-transfer coefficients, enhancement-factor parameters, packing area, liquid holdup, or hydraulic correlations. Its reaction coefficients describe equilibrium constants only.

**Unknown:** The current fitted equilibrium model's impact on absorber capture, temperature bulge, solvent circulation, pressure drop, and reboiler duty has not been quantified because the fitted package has not yet been connected and solved as the full column property model.

## Recommended staged path

1. **Property contract test.** Before constructing a column, build representative liquid states spanning the intended absorber and stripper ranges and verify fugacity, Henry's constant, true-species fractions, enthalpy, density, viscosity, surface tension, and diffusivities, including units and finite values.
2. **Resolve the shared eNRTL/enthalpy path.** Decide whether the fork's derivative bug should be fixed in `idaes-pse` or whether a project-specific enthalpy implementation is appropriate. Re-enable heat-of-absorption data only after the result is validated against the source data.
3. **Align the absorber interface.** Use one canonical species naming scheme, add the liquid transport methods required by `MEAColumn`, and provide a separately validated vapor package. Keep the equilibrium fit parameters separate from transport and kinetic parameters.
4. **Regression redesign.** Replace the arbitrary 10,000 speciation multiplier with documented measurement-error weights or a sensitivity study. Screen eNRTL interaction pairs, use source/temperature holdouts, and retain only parameter combinations that are identifiable over the intended operating envelope.
5. **Column bring-up.** Run a one- or two-element isothermal absorber first, then enable energy balances, then add the full discretization and realistic mass-transfer/kinetic correlations. Check solver termination and physical bounds at every stage.
6. **Process uncertainty.** Propagate correlated uncertainty in predicted equilibrium properties through capture, rich/lean loading, temperature profile, solvent circulation, and regeneration duty. Use parameter covariance only as a provisional local approximation until profile/bootstrap/Bayesian checks are available.

## External model context

The IDAES MEA property documentation lists liquid enthalpy/heat-of-absorption, Henry's constant, viscosity, thermal conductivity, vapor and liquid diffusivity, surface tension, and reaction-rate properties as part of the MEA property-method set: [IDAES liquid-phase property methods](https://idaes-pse.readthedocs.io/en/1.13.1/reference_guides/model_libraries/power_generation/carbon_capture/mea_solvent_system/properties/liquid_prop.html).

The official IDAES MEA packed-column implementation directly uses liquid surface tension, liquid viscosity, and liquid diffusivity in its mass-transfer equations: [IDAES `MEAsolvent_column.py`](https://github.com/IDAES/idaes-pse/blob/main/idaes/models_extra/column_models/MEAsolvent_column.py).

The broader modeling literature also treats VLE, enthalpy, solution chemistry, and reaction-model consistency together, and propagates thermodynamic parameter uncertainty through absorber/stripper process models: [Thermodynamic modeling and uncertainty quantification of CO2-loaded aqueous MEA solutions](https://doi.org/10.1016/j.ces.2017.04.049). Equilibrium-data uncertainty matters directly because equilibrium sets the mass-transfer driving force: [Experimentally based evaluation of accuracy of absorption equilibrium measurements](https://doi.org/10.1016/j.egypro.2013.05.176).

## Bottom line

Use the current reaction-equilibrium fit as a provisional thermodynamic initializer and for equilibrium/speciation comparisons. Do not yet use the raw coefficient uncertainties or the eNRTL-on smoke-fit values as process-design uncertainty. The next technically meaningful milestone is a property-contract and enthalpy validation, followed by a small isothermal rate-based absorber test; only then should process-level uncertainty propagation begin.
