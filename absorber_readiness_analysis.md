# eNRTL joint-regression and absorber-readiness analysis

Date: 2026-09-01  
Branch: `codex/enrtl-joint-caloric-regression`

## Result

**verified:** Replacing the correlated `[1, 1/T, ln(T)]` coefficient basis with the exact physical coordinates `ln(K_ref)`, `Delta H_ref`, and constant `Delta Cp` at `T_ref = 353.15 K` reduces the fitted reduced-Hessian condition number from `5.03e5` to `32.96`. The new Hessian eigenvalues span `36.36` to `1198.64`; the numerical Hessian defect is resolved.

**verified:** The accepted six-parameter joint fit uses the 240 VLE points, 172 speciation residuals, and 66 Kim-Svendsen calorimetry points at 40 and 80 C. The 20 points at 120 C are retained only as a temperature holdout. The fit weight on each squared heat residual is `0.002`, selected by a bounded objective tradeoff rather than copied from a larger published parameterization.

**verified:** The joint fit improves experimental-pressure MAPE from `39.28%` to `37.95%` and heat MAE from `9.59` to `7.26 kJ/mol CO2` relative to the equilibrium-only fit. The equilibrium objective rises only `3.3%`, from `59.66` to `61.66`.

**verified:** The untouched 120 C data fail: MAE is `27.96 kJ/mol CO2`, RMSE is `31.83 kJ/mol CO2`, and bias is `-20.05 kJ/mol CO2`. The present constant-Delta-Cp, reaction-only parameterization should not be used for non-isothermal absorber or stripper energy predictions near 120 C.

**conclusion:** The physical coordinate change fixes the ill-conditioned Hessian without adding parameters. It does not make the available data statistically sufficient. Reaction heat capacities remain weak, and the high-temperature holdout shows model-form or missing-data error larger than the fitted-data residuals.

## Thermodynamic formulation

The reaction equilibrium constants now use

```text
ln K(T) = ln K_ref
        + Delta H_ref/(R*T_ref)*(1 - T_ref/T)
        + Delta Cp/R*(ln(T/T_ref) + T_ref/T - 1)

Delta H(T) = Delta H_ref + Delta Cp*(T - T_ref)
```

**verified:** This is an exact coordinate transformation of the former constant, inverse-temperature, and logarithmic-temperature expression. Across 313.15-393.15 K, the transformed and former equations agree to `1.6e-14` in `ln K` when initialized from the same coefficients.

**verified:** Over the fitted temperatures, the normalized condition number of the former basis is `5812.5`; the normalized physical basis condition number is `2.84`. The effect-scaled physical basis condition number is `4.16`.

**reference-backed:** The constant-Delta-Cp form follows the standard Gibbs/enthalpy/heat-capacity relation used for the reaction standard states in the modified eNRTL formulation. It gives each fitted coordinate a direct thermodynamic interpretation.

**verified:** The local activity-equilibrium constraints remain dimensionless and use true-species activities with the aqueous-infinite-dilution reference state. The two reduced reactions conserve charge, MEA, CO2, and water:

```text
2 MEA + CO2 <-> MEAH+ + MEACOO-
MEA + CO2 + H2O <-> MEAH+ + HCO3-
```

**inference:** This reduced chemistry is appropriate over the fitted absorber-loading range. It is not established for every very-low-loading or hot stripper state.

## Calorimetry data correction

**verified:** The previous local Kim-Svendsen CSV merged two independent experimental runs at each temperature and globally sorted them. Near-equal loadings from different runs were then treated as consecutive doses, creating artificial increments such as `0.124` to `0.12401`.

**verified:** The corrected CSV transcribes the six source series separately: 19 and 18 points at 40 C, 15 and 14 at 80 C, and 10 and 10 at 120 C, for 86 points total. Heat differences are now formed only within an experimental run and begin from loading `0.003`.

**reference-backed:** Kim and Svendsen report semi-differential heats integrated over individual CO2 doses, generally spanning about `0.03-0.07 mol CO2/mol amine`, and report an experimental uncertainty of `+/-2.2%`.

**reference-backed:** Akula et al. fit the 66 Kim-Svendsen points at 40 and 80 C and exclude the 120 C series because prior literature identified thermodynamic inconsistency at and above 373 K. Their objective heat weight of `0.2` belongs to a larger ten-parameter standard-state and ion-pair regression, not this six-parameter reaction-only fit.

### Independent Kim et al. (2014) validation

**verified:** Kim et al. Table A1-1 on p. 1455 contains 27 semi-differential observations for 30 wt% MEA: 11 at 40 C, 9 at 80 C, and 7 at 120 C. All 27 are retained with `fit=0`; none influence the fitted parameters.

**verified:** The current model gives MAE values of `3.25`, `3.49`, and `12.26 kJ/mol CO2` at 40, 80, and 120 C, respectively. The corresponding biases are `+1.15`, `-2.22`, and `-12.26 kJ/mol CO2`. Every 120 C prediction is below its observation.

**verified:** Kim et al. report temperature-sensor accuracy, pressure-transducer accuracy, calorimeter sensitivity, and agreement within 5% between loading from liquid analysis and flowmeter readings. They do not report a direct percentage uncertainty for heat of absorption, so the retained 2014 rows have no invented heat-error bars.

**inference:** Agreement at 40 and 80 C rejects a general failure of the semi-differential enthalpy expression. The one-sided 120 C residual identifies high-temperature caloric temperature dependence as the remaining model or data disagreement.

**verified:** Kim et al. state on p. 1448 that their 80 and 120 C heats are approximately 5% and 10% lower than the earlier 2007 values and found no reasonable explanation. The 2014 result therefore narrows the discrepancy but does not establish that either high-temperature source is unbiased.

## Weight selection

| Heat weight | Pressure MAPE | Heat MAE (fit) | Hessian condition | Decision |
|---:|---:|---:|---:|---|
| `0` | `39.28%` | `9.59 kJ/mol` | `36.61` | equilibrium-only baseline |
| `0.002` | `37.95%` | `7.26 kJ/mol` | `32.96` | accepted |
| `0.01` | `45.85%` | `5.68 kJ/mol` | `67.13` | pressure degradation |
| `0.2` | `137.13%` | `3.67 kJ/mol` | `2138.62` | incompatible with reduced parameter set |

**verified:** At weight `0.002`, the heat MAE is `4.75 kJ/mol` at 40 C and `10.46 kJ/mol` at 80 C. The 80 C bias is `-10.17 kJ/mol`, showing that the accepted compromise still has a systematic temperature trend.

**inference:** Increasing the heat weight cannot repair that trend with the current parameter set; it rotates the regression away from the VLE data and drives the reaction heat capacities to compensate.

## Parameters and curvature

| Reaction | Coordinate | Value | Local curvature scale | Relative scale |
|---|---|---:|---:|---:|
| bicarbonate | `ln K_ref` | `5.2222` | `0.1379` | `2.64%` |
| bicarbonate | `Delta H_ref` | `-67.536 kJ/mol` | `3.774 kJ/mol` | `5.59%` |
| bicarbonate | `Delta Cp` | `-167.82 J/mol/K` | `116.97 J/mol/K` | `69.7%` |
| carbamate | `ln K_ref` | `11.6819` | `0.1119` | `0.96%` |
| carbamate | `Delta H_ref` | `-102.786 kJ/mol` | `3.159 kJ/mol` | `3.07%` |
| carbamate | `Delta Cp` | `-241.98 J/mol/K` | `120.96 J/mol/K` | `50.0%` |

**verified:** The largest absolute Hessian-derived correlation is `0.759` between bicarbonate `Delta H_ref` and `Delta Cp`; the corresponding carbamate correlation is `0.646`.

**verified:** The output columns are now named `Curvature_Scale` and `Relative_Curvature_Scale`. They are not called uncertainty because the VLE, speciation, and calorimetry residuals do not have a common calibrated measurement covariance, and the objective contains regularization.

**unknown:** Frequentist confidence intervals or a posterior parameter covariance cannot be justified from the current objective. A source-level error model and additional independent calorimetry are required before making those claims.

## eNRTL interactions and absorber use

**reference-backed:** Akula et al. simultaneously regress ion standard-state formation properties, heat-capacity terms, and selected water-ion-pair interactions. Their reported parameter table contains severe correlations and very large variances for some heat-capacity and reverse-interaction terms.

**verified:** Earlier attempts to add 12 local-interaction variables to the available objective produced a non-positive-definite reduced Hessian. Applying the published fixed interactions without the accompanying standard-state regression produced a much worse objective (`1333.36`).

**conclusion:** Adding eNRTL interaction parameters now would add ambiguity rather than reduce it. The six physical reaction coordinates are the smallest identifiable regression supported by the retained local data.

**verified:** The local package still lacks the liquid diffusivity, viscosity, surface tension, thermal conductivity, and vapor-phase interfaces needed by the rate-based absorber. It also uses ion names different from the sibling column package. These interface gaps are separate from the equilibrium Hessian fix.

## Recommended next scientific step

1. Profile the two existing reaction `Delta Cp` coordinates against the 2014 validation residuals while refitting the other coordinates to the original objective.
2. Retain a shift only if it improves the independent 120 C comparison without materially degrading VLE, speciation, or the independent 40 and 80 C comparisons.
3. Add a standard-state heat-capacity contribution only if neither existing reaction coordinate can satisfy that comparison and the new coordinate remains identifiable.
4. Propagate the resulting correlated prediction range through a small isothermal column case before enabling the energy balance in a full absorber.

Do not fit the 120 C Kim-Svendsen holdout merely to lower its error. Its failure is useful evidence that the current caloric temperature dependence is not transferable.

## Reproduction

```bash
env MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 .venv/bin/python Fitting_Routine.py
.venv/bin/python Calorimetry_Validation.py generate
.venv/bin/python Calorimetry_Validation.py render
```

The retained comparison values are in `data/Plots/Calorimetry_Validation.csv`; the model-versus-observation figure is `data/Plots/Calorimetry_Validation.png`. The fit also writes `data/Parameters/Parameter_Correlation.csv`.

## Primary references

- Kim and Svendsen, “Heat of Absorption of Carbon Dioxide (CO2) in Monoethanolamine (MEA) and 2-(Aminoethyl)ethanolamine,” *Industrial & Engineering Chemistry Research* 46, 5803-5809, 2007, [doi:10.1021/ie0616489](https://doi.org/10.1021/ie0616489).
- Kim, Hoff, and Mejdell, “Heat of Absorption of CO2 with Aqueous Solutions of MEA: New Experimental Data,” *Energy Procedia* 63, 1446-1455, 2014, [doi:10.1016/j.egypro.2014.11.154](https://doi.org/10.1016/j.egypro.2014.11.154).
- Akula et al., “A Modified Electrolyte Non-Random Two-Liquid Model with Analytical Expression for Excess Enthalpy: Application to the MEA-H2O-CO2 System,” *AIChE Journal* 69(1), 2023, [doi:10.1002/aic.17935](https://doi.org/10.1002/aic.17935).
