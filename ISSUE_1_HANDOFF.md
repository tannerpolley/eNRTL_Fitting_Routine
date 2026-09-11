# Issue #1 implementation checkpoint

This is incomplete work, stopped at the parent's placement handoff on 2026-09-11.
No push, PR, corrected fit, or immutable export has been completed. Existing
fit/calorimetry/profile/absorber tables remain historical and must not be
presented as corrected results.

Base: `af3fade`. Branch: `codex/issue-1-corrected-mea-packet`.
Parent approved a ceiling of 650 gross executable additions against this base;
tests, documentation, and generated tables are accounted separately.

## Completed changes

- `ENRTLDirectTemperatureDerivative._dA_DH_dT` includes explicit `3/T`;
  the excess-enthalpy owner uses it. The old expression remains only in the test.
- Fixed-true-composition central differences at 313.15, 353.15, and 393.15 K,
  with steps 1, 0.1, 0.01, 0.001, and 0.0001 K, pass the corrected expression
  and reject the old expression. Values are retained in
  `data/Plots/Debye_Huckel_Derivative.csv`.
- `get_prop_dict` no longer consumes the global reaction combination dictionary,
  allowing repeated model construction in one process.
- `Profile_DeltaCp.equilibrium_metrics` deterministically pools residuals and
  records per-source counts, units, and explicit fraction/percent MAPE.
  Its test uses unequal source sizes and reversed source order.
- Dataset iteration is sorted. Plot-only aggregate MAPE is weighted by row count
  and labeled as including Xu; it is distinct from direct calibration metrics.
- The fitting script now retains direct residuals and fit metrics. Hessian
  analysis now retains the newly evaluated matrix and eigenvalues.
- Profiles read curvature scales from the newly generated parameter table;
  they retain pooled/source-separated metrics and parameters for each point.
- A fixed-bicarbonate run writes separate historical-coordinate sensitivity
  filenames. No new candidate selection rule was introduced.
- Kim 2014 is labeled selection-exposed comparison.

## Validation and interrupted execution

Two focused unittest checks passed in approximately 0.2 seconds (excluding
imports), using the canonical fitting checkout's existing Python environment.
The installed IDAES distribution imports successfully and its direct URL
identifies `dallan-keylogic/idaes-pse` commit
`3cdf9b5299b8453fc29274f55d060c1b30a8294b`, requested revision `eNRTL`.
Python is 3.12.13. Solvers are in the existing IDAES bin directory.
No dependencies or downstream repositories were modified.

The full `Fitting_Routine.py` run was started with `MPLBACKEND=Agg` and
`PYTHONDONTWRITEBYTECODE=1`. At the handoff it had reported `DOF: 0` and was
still initializing/building; it had not produced a solver result, corrected
parameter table, or fresh curvature. The parent-requested stop sent SIGTERM
only after checking the exact process command and worktree directory. This is
an interrupted execution, not numerical nonconvergence or a model failure.

## Remaining work, in order

1. Review this checkpoint and run the frozen full baseline fit to completion.
   Preserve the objective weights, six reaction parameters, excluded Xu and
   Bottinger datasets, loading filter, starts, regularization, solver bounds,
   tolerances, and iteration budget. Diagnose any failure before changing policy.
2. Regenerate calorimetry, fitted plots, fresh Hessian/correlation/scales,
   ten-point Delta-Cp profiles, and (if retained) the historical-coordinate
   sensitivity. Run the profile only after the corrected fit overwrites the
   old parameter table. Preserve unfavorable 120 C outcomes.
3. Check actual observation counts and direct metric aggregation; test the
   production enthalpy routing in addition to the current A-derivative check,
   plus repeatable configuration and remaining acceptance criteria.
4. Replace or clearly supersede stale conclusions in
   `absorber_readiness_analysis.md`. The canonical uncommitted Cp additions
   were inspected and deliberately not copied: apparent ideal-mixture Cp does
   not establish full reacting-solution Cp. Total solution Cp stays unavailable.
5. Do not run the absorber adapter or change either downstream repository.
   Retain historical process evidence only with exact byte/commit identity and
   its incomplete-runtime/model-to-model limitations.
6. Implement the compact immutable export and validator: full thermodynamic
   equations/conventions, species/reaction ordering and direction, all numerical
   parameters and units, reference temperatures/states, Henry/pressure and
   electrostatic/local terms, observation roles/domains, dependency/solver
   identity, results, limitations and unavailable fields. Include exact-byte
   SHA-256 inventory and external manifest digest; exclude absolute paths and
   oversized files. No export generator or schema has yet been written.
7. Complete focused packet/acceptance tests, size check, cleanup audit,
   independent acceptance, and final change accounting. Commit, push, and open
   a draft PR closing issue #1 only under the continuing parent's authority;
   do not merge.

Reproduction commands, with `PYTHON` set to the compatible existing environment:

```sh
PYTHONDONTWRITEBYTECODE=1 "$PYTHON" -m unittest discover -s tests -v
MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 "$PYTHON" Fitting_Routine.py
MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 "$PYTHON" Calorimetry_Validation.py generate
MPLBACKEND=Agg "$PYTHON" Calorimetry_Validation.py render
MPLBACKEND=Agg "$PYTHON" Profile_DeltaCp.py run
MPLBACKEND=Agg "$PYTHON" Profile_DeltaCp.py render
MPLBACKEND=Agg "$PYTHON" Profile_DeltaCp.py candidate
```

The canonical checkout was inspected read-only. Its uncommitted files remain
untouched. A replacement task should reuse this branch/checkpoint rather than
reconstructing changes from the canonical dirty checkout.
