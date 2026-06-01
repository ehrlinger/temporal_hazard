# TemporalHazard 1.1.0 (development version)

## Bug fixes

* **4-phase CoE fixmu-phase selection** (Phase 7d).
  `.hzr_select_fixmu_phase()` used `which.max()` over raw per-phase cumhaz
  at the starting theta.  G3 late phases with typical shape parameters have
  unnormalized cumhaz orders of magnitude larger than other phases, causing
  CoE to pin the G3 `log_mu` away from its true near-zero MLE.  Fixed by
  excluding phases whose cumhaz contribution exceeds 10× the median before
  selecting (falls back to `which.max` when all phases are outliers).  On the
  4-phase CABGKUL fit the CoE vs no-CoE LL gap closes from 6.9 to < 0.1
  units.  Six new tests cover the 4-phase code path in
  `test-conservation-of-events.R`.
* **`time_lower` dual-use bug in Weibull and multiphase likelihoods.**
  When `time_lower` was supplied for a mixed interval-censored + right-censored
  dataset, the Weibull LL interpreted `time_lower` as the counting-process
  *entry time* for right-censored rows, computing H(stop) − H(start) = 0 and
  silently zeroing those rows' likelihood contribution.  Fixed in
  `likelihood-weibull.R` (4 sites: LL, gradient, L-BFGS-B internal LL/gradient)
  and `likelihood-multiphase.R`: `start_vec` is now set from `time_lower` only
  for genuine epoch rows (`status %in% c(0L, 1L)` and `time_lower < time`).
  Two regression tests added to `test-interval-censoring-weibull.R`.

* **`hzr_decompos()` Case 3 corrected and `nu = 0, m >= 0` now fails loud**
  (Phase 7d).  Two issues in the early-phase (G1) sign dispatch:
    - **Case 3 (`m > 0, nu < 0`, "bounded cumulative") carried a spurious
      factor of `m`.**  Its `rho` used a bare `(2^m - 1)^nu` instead of the
      `((2^m - 1)/m)^nu` form used by Case 1, leaving an `m` factor on the
      `bt^(-1/nu)` term.  The CDF diverged from the C HAZARD G1 evaluator
      (`g1flag = 5`) by up to ~0.2 and was discontinuous with its `m -> 0`
      limit (Case 3L).  Adding the `/m` divisor makes the `m` factors cancel,
      reproducing the C evaluator exactly and restoring continuity (verified
      against `src/common/hzd_ln_G1_and_SG1.c`).  No shipped phase uses
      Case 3, so fitted models are unaffected; the synthetic 3-phase golden
      fixture was regenerated because its free-shape optimizer path crosses
      Case 3 territory.
    - **`nu = 0` with `m >= 0`** fell through every dispatch branch, leaving
      the CDF unassigned and raising the cryptic `object 'G' not found`.  The
      `nu -> 0` limit is defined only for `m < 0`; for `m >= 0` it is
      degenerate.  The function now raises a clear, explanatory error.
  New `test-decompos-boundary.R` locks in continuity of all limiting branches
  (Case 1 -> 1L, 2 -> 1L, 2 -> 2L, 3 -> 3L), Case 3 <-> C `g1flag=5` parity,
  `g = dG/dt` internal consistency, CDF sanity, and stability at extreme
  `t_half`.

## New features

* `vignette("fitting-hazard-models")` gains an **Interval and left censoring**
  section covering: status coding reference (`-1`/`0`/`1`/`2`), a cardiac
  clinic-visit simulation with right- and interval-censored observations,
  the direct `time_lower`/`time_upper` API, and a comparison showing the
  interval-censored fit recovering `nu` close to 1.0 (true value) while the
  naive exact-at-upper fit incurs a shape bias of ~+0.45.  Includes a callout note on the correct
  use of `time_lower = 0` for right-censored rows.
* `vignette("fitting-hazard-models")` gains a **Convergence troubleshooting**
  section covering: reading the KM cumulative hazard for Weibull starting
  values (log-log plot), when to fix shape parameters vs. estimate freely,
  diagnosing overparameterization via near-zero phase scales and `NA` from
  `vcov()`, and `control` options (`n_starts`, `maxit`).

---

# TemporalHazard 1.0.3

## Bug fixes / CRAN compliance

* `hzr_bootstrap()` no longer touches `.GlobalEnv` directly. The 1.0.2
  `oldseed`/`on.exit()`/`assign(".Random.seed", ...)` save-restore wrapper
  added in 1.0.2 violated CRAN policy on writing to `.GlobalEnv` and has
  been removed. When `seed` is supplied the function simply calls
  `set.seed(seed)` (the documented R API for seeded reproducibility); the
  `@param seed` documentation now notes that the caller's RNG state is not
  restored on exit. With `seed = NULL` (the default) the function does
  not call `set.seed()` at entry, so it starts from the caller's current
  RNG state; the bootstrap still consumes random numbers and advances
  that state in the usual way.

# TemporalHazard 1.0.2

## Bug fixes / CRAN compliance

* The golden-fixture generators (`.hzr_create_*_golden_fixture()`,
  previously `R/golden_fixtures.R`) have been moved out of the package to
  `data-raw/golden_fixtures.R`. They are maintainer-only helpers for
  regenerating the bundled `inst/fixtures/*.rds` reference outputs and are
  not part of the installed package, so they are no longer shipped, checked,
  or user-reachable. This resolves the home-filespace concern at its root:
  the earlier fallback resolved to `system.file("fixtures", ...)` — i.e. the
  installed package directory — whenever the package was installed, so the
  1.0.1 "falls back to `tempdir()`" fix did not actually prevent writing to
  the user library. The bundled `.rds` fixtures still ship and the parity
  tests still read them via `system.file()`.
* `.hzr_generate_golden_fixture()` (the C-binary reference writer in
  `R/parity-helpers.R`, which shares a file with test-time helpers and so
  was kept in the package) now takes a required `output_dir` argument with
  no default path.
* Removed the remaining hardcoded `seed = 42` literals from the relocated
  generators; recorded fixture metadata reflects the actual `seed` argument
  passed (`NULL` by default, so no seed is set inside the function).
* `hzr_bootstrap()` no longer leaves the caller's random-number stream
  altered when `seed` is supplied: the global `.Random.seed` is saved before
  `set.seed()` and restored via `on.exit()`, matching the fixture generators.
  Bootstrap reproducibility under a given `seed` is unchanged.

# TemporalHazard 1.0.1

## Bug fixes / CRAN compliance

* Added `\value` documentation to all exported functions that were missing it:
  `hazard()`, `coef.hazard()`, `vcov.hazard()`, `print.hzr_calibrate()`,
  `print.hzr_deciles()`, `print.hzr_gof()`, and `print.hzr_kaplan()`.
* Internal fixture generators (`R/golden_fixtures.R`) no longer set a specific
  seed unconditionally. Generators now accept an optional `seed` argument;
  when provided, the global RNG state is saved and restored via `on.exit()`.
* Default `output_dir` for fixture generators falls back to `tempdir()` instead
  of the package source directory, keeping the home filespace unmodified.

# TemporalHazard 0.9.8

## New features

* **Delta-method confidence limits on `predict.hazard()`** — Phase 4g of
  the development plan lands. Two new arguments: `se.fit = FALSE` and
  `level = 0.95`. When `se.fit = TRUE`, the return value becomes a
  data frame with columns `fit`, `se.fit`, `lower`, `upper`.
  - **Weibull and multiphase use closed-form Jacobians**
    (`dH/dtheta`, `dexp(eta)/dtheta`, `deta/dtheta`); exponential /
    log-logistic / log-normal fall back to `numDeriv::jacobian` on a
    per-call cumhaz closure.
  - **Transforms match SAS HAZARD** (`hzp_calc_haz_CL.c` /
    `hzp_calc_srv_CL.c`): `hazard` and `cumulative_hazard` use
    log-scale CLs; `survival` uses log(-log S) CLs (equivalent to
    log-cumhaz) so 0 <= lower <= upper <= 1; `linear_predictor` is
    symmetric on the natural scale.
  - **Fixed-shape / CoE multiphase fits produce meaningful CLs** — the
    delta-method sandwich is restricted to the free-parameter submatrix
    of `vcov`, treating fixed parameters as known-with-zero-variance.
  - Backward compatible: `se.fit = FALSE` (default) preserves the
    pre-0.9.8 scalar-vector / decompose-data-frame return shape.

# TemporalHazard 0.9.7

## New features

* **Counting-process / repeating-events likelihood wired up** — Phase 4f
  of the development plan lands. `Surv(start, stop, event)` with any
  `start > 0` is now accepted. The Weibull and multiphase log-likelihoods
  apply `H(stop) - H(start)` to event and right-censored terms; the
  trivial `start = 0` case degenerates to `H(stop)` and recovers the
  plain-Surv fit exactly. Splitting each row into contiguous epochs
  preserves both the log-likelihood and the MLE to optimizer tolerance
  (split-invariance).
* **Weibull + multiphase analytic gradients handle H(start).** The
  closed-form Weibull score adds a `-d H(start)/d theta` term per row
  (guarded at `start = 0`). The multiphase analytic gradient computes
  per-phase `Phi_j(start)` and its shape derivatives, then adds
  `+w_H_start * mu_j * dPhi_j(start)` to each parameter's score; G3
  phase derivatives at `start` use the same finite-difference machinery
  as at `stop`.
* **0.9.5 narrowing removed.** The `hazard()` guard that rejected
  counting-process `Surv(start, stop, event)` with any `start > 0` is
  gone.

# TemporalHazard 0.9.6

## New features

* **`weights` now supported for all distributions** — Phase 4e of the
  development plan lands. The exponential, log-logistic, and
  log-normal likelihoods and their analytic gradients now apply row
  weights to every censoring term (event, right-censored,
  left-censored, interval-censored). The 0.9.5 guard in `hazard()`
  that rejected `weights` for `dist %in% c("exponential",
  "loglogistic", "lognormal")` has been removed. Fits with integer
  weights reproduce the row-duplicated fit to optimizer tolerance
  across all five distributions.
* **Conservation of Events now honours weights.**
  `.hzr_conserve_events()` and `.hzr_select_fixmu_phase()` take an
  optional `weights` argument; the multiphase optimizer threads it
  through so per-phase cumulative hazards are summed on the same
  scale as the (weighted) observed event count. CoE no longer
  auto-disables when weights are non-uniform — the dimension
  reduction stays on and the MLE matches the full-dim path.

## Bug fixes

* **Multiphase analytic gradient now applies `weights`.**
  `.hzr_gradient_multiphase()` accepted neither `weights` nor its
  downstream equivalents: the per-row score weights `w_H` / `inv_h`
  were set to ±1 and the interval-censored finite-difference
  correction summed an unweighted LL. Weighted multiphase fits
  therefore optimised a weighted objective with an unweighted score;
  BFGS line search still converged near the correct MLE but the
  final gradient norm did not go to zero. All three paths now honour
  row weights, and the optimizer's `gradient_fn` wrapper (including
  the all-zero numeric fallback and the CoE wrapper) forwards
  `weights` consistently. Regression test covers weighted analytic
  vs numerical gradient parity. Surfaced by Copilot review on PR #18.

# TemporalHazard 0.9.5

## New features

* **Stepwise covariate selection** — `hzr_stepwise()` runs forward,
  backward, or two-way stepwise selection on an existing `hazard` fit
  using Wald p-values or AIC deltas as the entry / retention criterion.
  Phase-specific entry is supported for multiphase models: a covariate
  can enter one phase and not another. Defaults match SAS `PROC HAZARD`
  (`SLENTRY = 0.30`, `SLSTAY = 0.20`); AIC mode uses `ΔAIC < 0`
  uniformly. SAS-style `MOVE` oscillation guard freezes variables that
  enter + exit more than `max_move` times. Returns an object of class
  `c("hzr_stepwise", "hazard")` with a `$steps` selection trace, scope
  record, and elapsed timer. Implements the core algorithm from C
  HAZARD `stepw.c` / `backw.c`.

## Bug fixes

* **Multiphase convergence after weights/repeating-events merge** —
  restored multiphase optimization that regressed in 0.9.4: three
  interacting defects in the new `weights` threading (dup-arg
  collision in the multiphase / Weibull closures, positional-arg
  corruption in every distribution's gradient call) made every
  optimizer iteration error silently inside `tryCatch`. Diagnosed and
  fixed via commit 73b4657.
* **Weibull analytic gradient now applies `weights`** — both
  `.hzr_gradient_weibull()` and the `grad_internal` closure inside
  `.hzr_optim_weibull()` accepted `weights` as a formal but did not
  apply it to the score vector. The optimizer still converged via
  line search on the (weighted) log-likelihood, but the gradient
  direction was wrong and the final gradient norm did not go to
  zero. Both gradient paths now weight the event indicator and
  cumulative hazard building blocks. Fits with integer weights
  reproduce the equivalent row-duplicated fit to optimizer tolerance.

## Scope change

* `weights` is now only accepted for `dist = "weibull"` and
  `dist = "multiphase"`. The 0.9.4 NEWS claimed weights were
  threaded through all distribution-specific likelihoods; in fact the
  exponential, log-logistic, and log-normal single-distribution paths
  accepted the formal but never applied it, so the fit was silently
  unweighted. `hazard()` now raises an explicit error when `weights`
  is supplied with one of those distributions rather than returning
  an unweighted fit. Full support for the remaining single-dist paths
  is tracked in `inst/dev/DEVELOPMENT-PLAN.md` Phase 4e.
* **Conservation of Events is auto-disabled when weights are not
  all 1.** `.hzr_conserve_events()` receives the weighted event count
  as its target but sums per-phase cumulative hazards across rows
  *without* applying weights, so Turner's adjustment comes out on a
  mismatched scale. The multiphase optimizer now detects non-unit
  weights and skips the CoE dimension reduction, falling through to
  the (correctly weighted) full-dimensional path. Fits are still
  correct; they just don't benefit from the one-parameter
  analytical closed-form solve. Weighted CoE wire-up is tracked
  alongside the other weights completion work in
  `inst/dev/DEVELOPMENT-PLAN.md` Phase 4e.
* **Repeating-events / counting-process notation narrowed.**
  `Surv(start, stop, event)` with `start > 0` is no longer accepted
  by `hazard()`. The 0.9.4 NEWS claimed each epoch contributed
  `H(stop) - H(start)` to the likelihood, but downstream likelihoods
  only read `time_lower` for interval-censored rows (`status == 2`);
  counting-process rows (`status` in `{0, 1}`) were silently scored
  with `H(stop)` alone, so any fit with nonzero entry times was
  silently wrong. `hazard()` now raises an explicit error. The
  trivial case `Surv(0, t, d)` -- equivalent to `Surv(t, d)` --
  continues to work. Full wire-up of `H(stop) - H(start)` for all
  distribution paths is tracked in
  `inst/dev/DEVELOPMENT-PLAN.md` Phase 4f.

# TemporalHazard 0.9.4

## New features

* **Observation weights** — `weights` argument in `hazard()` applies Fisher
  weighting to the log-likelihood for `dist = "weibull"` and
  `dist = "multiphase"`. Each observation's contribution is multiplied
  by its weight, enabling severity-weighted event analyses. Implements
  the SAS `WEIGHT` statement. _The original 0.9.4 entry claimed
  coverage of all distribution paths; the 0.9.5 patch corrected the
  claim and fixed a gradient wire-up bug in the Weibull path._
* **Repeating events** — `Surv(start, stop, event)` start-stop notation
  is parsed. _The original 0.9.4 entry claimed each epoch contributed
  `H(stop) - H(start)` to the likelihood, but the downstream
  likelihoods never applied the lower bound for counting-process rows;
  the 0.9.5 patch narrowed the feature to the trivial `start = 0`
  case and added an explicit error for nonzero starts._

# TemporalHazard 0.9.3

## New features

* `hzr_deciles()` — Decile-of-risk calibration function comparing observed
  vs. expected event counts across risk groups with chi-square GOF testing.
  Implements the SAS `deciles.hazard.sas` macro workflow.
* `hzr_gof()` — Goodness-of-fit function comparing parametric predictions
  against nonparametric (Kaplan-Meier) estimates with observed vs. expected
  event counting. Implements the SAS `hazplot.sas` macro workflow.
* `hzr_kaplan()` — Kaplan-Meier survival estimator with logit-transformed
  confidence limits that respect the [0, 1] boundary, interval hazard rate,
  density, and restricted mean survival time (life integral). Implements the
  SAS `kaplan.sas` macro output structure.
* `hzr_calibrate()` — Variable calibration function for assessing functional
  form before model entry. Groups a continuous covariate into quantile bins
  and applies logit, Gompertz, or Cox link transforms. Supports
  stratification via the `by` parameter. Implements the SAS `logit.sas` and
  `logitgr.sas` macros.
* `hzr_nelson()` — Wayne Nelson cumulative hazard estimator with lognormal
  confidence limits. Supports weighted events for severity-adjusted repeated
  event analyses. Implements the SAS `nelsonl.sas` macro.
* `hzr_bootstrap()` — Bootstrap resampling for hazard model coefficients with
  bagging support (fractional sampling). Returns per-replicate estimates and
  summary statistics (mean, SD, percentile CI). Implements the SAS
  `bootstrap.hazard.sas` macro workflow.
* `hzr_competing_risks()` — Competing risks cumulative incidence using the
  Aalen-Johansen estimator with Greenwood variance. Handles any number of
  competing event types. Implements the SAS `markov.sas` macro.
* **Conservation of Events (CoE)** — Turner's theorem is now integrated into
  the multiphase optimizer. One phase's log_mu scaling parameter is solved
  analytically at each iteration, reducing the optimization dimension by 1
  and improving numerical stability and convergence. Enabled by default;
  disable with `control = list(conserve = FALSE)`. Implements the core
  algorithm from C HAZARD `setcoe.c` / `consrv.c`.
* New vignette: "Complete Clinical Analysis Walkthrough" — end-to-end
  workflow from Kaplan-Meier baseline through validated multivariable model,
  mirroring the SAS HAZARD analytical sequence.

## Improvements

* Multi-start optimizer now respects user-set RNG seeds for reproducibility
  (removed `set.seed(NULL)` that was actively breaking determinism).
* Vignette metadata normalized to YAML `vignette:` key across all 8 files.
* `fit` parameter documentation corrected to state default is FALSE.
* README now includes key capabilities table and development plan link.

# TemporalHazard 0.9.1

## New features

* G3 late-phase decomposition (`hzr_phase("g3", ...)`) now fully integrated
  into the multiphase optimizer, Hessian, and prediction pipeline.
* `fixed = "shapes"` parameter in `hzr_phase()` allows fixing shape parameters
  during estimation (matching C/SAS HAZARD workflow of estimating only log-mu
  scale parameters).

## Bug fixes

* `summary.hazard()` now correctly reports standard errors when some
  parameters are fixed. Previously, `anyNA(vcov)` rejected the entire
  variance-covariance matrix when fixed parameters had NA entries.
* `print.summary.hazard()` coefficient table now shows the correct label
  for G3 phases (was printing empty parentheses).
* `print.summary.hazard()` phase listing now uses the phase name in
  CDF labels (e.g., "cdf (late risk)") instead of hardcoded "early risk".
* SAS missing value markers (`.`) in CSV datasets are now handled via
  `na.strings = c("NA", ".")` in `data-raw/make_data.R`, preventing
  numeric columns from being read as character.

## Documentation

* Seven Quarto vignettes: getting-started, fitting-hazard-models,
  prediction-visualization, inference-diagnostics, mathematical-foundations,
  package-architecture, and sas-to-r-migration.
* Roxygen examples now include both single-phase and multiphase models.
* README switched to self-contained CABGKUL examples with G3 late phase.
* Dataset axis labels corrected to "Months" (not "Years").

## Infrastructure

* CI workflows updated to use `roxygen2::load_pkgload` for lazy data
  compatibility.
* Added lintr CI workflow with `.lintr` configuration.
* pkgdown action bumped to `peaceiris/actions-gh-pages@v4`.
* Added `use-public-rspm: true` to all CI workflows.
* Added `lintr` to Suggests.

# TemporalHazard 0.9.0

## New features

* Multiphase engine: N-phase additive cumulative hazard models via
 `dist = "multiphase"` with `hzr_phase()` specification.
* `hzr_decompos()` parametric family implementing the three-parameter
  temporal decomposition of Blackstone, Naftel, and Turner (1986).
* Multi-start optimizer with Hessian-based variance-covariance estimation.
* C binary parity tests against the KUL CABG reference dataset.
* Five clinical reference datasets: `avc`, `cabgkul`, `omc`, `tga`, `valves`.

# TemporalHazard 0.1.0

## New features

* Single-phase engine: Weibull, exponential, log-logistic, and log-normal
  distributions with formula interface.
* `hazard()` API with `predict()`, `summary()`, `coef()`, `vcov()` S3 methods.
* Golden fixture regression testing system.
* Numerically stable helper primitives (`hzr_log1pexp`, `hzr_log1mexp`,
  `hzr_clamp_prob`).

# TemporalHazard 0.0.0.9000

* Initial package scaffold.
* Added numerically stable helper primitives.
* Added baseline unit tests and CI workflow.
