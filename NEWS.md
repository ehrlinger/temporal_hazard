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

# TemporalHazard 0.9.4

## New features

* **Observation weights** — `weights` argument in `hazard()` applies Fisher
  weighting to the log-likelihood. Each observation's contribution is
  multiplied by its weight, enabling severity-weighted event analyses.
  Threaded through all distribution-specific likelihood and gradient
  functions. Implements the SAS `WEIGHT` statement.
* **Repeating events** — `Surv(start, stop, event)` start-stop notation
  is now supported for epoch-decomposed longitudinal data. Each epoch
  contributes `H(stop) - H(start)` to the likelihood, matching the SAS
  `LCENSOR`/`STIME` mechanism.

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
