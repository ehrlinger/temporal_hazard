# CRAN submission comments -- TemporalHazard 1.1.0

## Summary

This is an update to TemporalHazard 1.0.3 (accepted 2026-05-29).
No reviewer feedback to address. New version adds analytic Hessians for all
five distribution families, full-information variance for
Conservation-of-Events fits, and a suite of SAS HAZPRED parity fixtures.

## Changes since 1.0.3

* **Analytic Hessian for all five families (Phase 7c, Layer 2).** The
  exponential, Weibull, log-logistic, log-normal, and multiphase
  distributions now compute post-fit Hessians in closed form, replacing
  the numerical Richardson approximation used previously. Each analytic
  Hessian is validated against `numDeriv` at a non-trivial parameter
  point; left/interval-censored fits fall back to the numerical Hessian.
  The multiphase Hessian covers the three-term NLL assembly including the
  Conservation-of-Events full-information path.

* **Full-information variance for Conservation-of-Events fits.** CoE
  constrains one phase's `log_mu` during optimization but previously also
  dropped it from the uncertainty, producing an `NA` standard error for
  the conserved phase. The conserved-phase SE is now recovered by
  recomputing the unconstrained Hessian over the full free set at the
  optimum. On `hz.death.AVC` every SE now matches the SAS HAZARD
  reference.

* **`predict.hazard()` improvements.** Added `type = "hazard"` for
  multiphase (instantaneous additive hazard), `conf.type` for survival
  confidence-limit transform (log-log default, logit for SAS HAZPRED
  parity), and `decompose = TRUE` + `se.fit = TRUE` for per-phase
  delta-method confidence limits on multiphase cumulative hazard.

* **`hzr_deciles()` now matches the SAS `deciles.hazard` macro exactly.**
  Subjects are ranked into equal-sized risk groups by predicted survival
  at the horizon; expected counts use the sum of predicted cumulative
  hazard at each subject's own follow-up time.

* **Hessian inversion hardening.** Post-fit variance-covariance now
  symmetrizes the Hessian, checks conditioning via `rcond`, applies a
  Cholesky-with-`solve()` fallback, and guards against non-positive
  variances. Ill-conditioned fits carry `rcond` and `pd` diagnostics
  surfaced in `summary()`.

* **SAS HAZPRED parity fixtures.** Group A fixture suite covers patient-
  specific predictions, stratified calibration, Kaplan-Meier / Nelson-
  Aalen life tables, and decile-of-risk tables for the `hz.death.AVC`
  reference model; all verified to print-precision or better.

* **Formula preprocessor hardened.** The multiphase formula preprocessor
  now uses a parse-tree walk to detect phase-scoped calls, replacing a
  `deparse`+regex approach that produced false positives when a phase name
  coincided with a base-R function name.

* **`hzr_bootstrap()` robustness.** Resamples `object$data$weights` and
  `object$data$frame` directly rather than re-evaluating the call in the
  caller frame, so bootstraps work correctly when the original data or
  weight objects are out of scope.

## Test environments

* **Local:** R 4.6.0 on macOS (aarch64-apple-darwin23).
  `R CMD check --as-cran` (with PDF manual) returns 0 errors, 0 warnings,
  0 notes (on the release-version tarball; the `.9000` dev tarball shows
  an expected "Version contains large components" NOTE that disappears
  once `.9000` is stripped).
* **GitHub Actions matrix:** ubuntu-latest (R-devel / R-release /
  R-oldrel-1), macos-latest (R-release), windows-latest (R-release).
  All checks passed (0 errors / 0 warnings / 0 notes).
* **Reverse-dependency check:** `tools::package_dependencies(reverse = TRUE)`
  returns 0.

## NOTE disposition

`R CMD check --as-cran` on the release tarball is clean (0/0/0).
CRAN's incoming-feasibility / aspell check may report:

* **Possibly misspelled words in DESCRIPTION:** `Naftel`, `Rajeswaran`
  (cited authors with `<doi:...>` references), `UAB` (acronym), `et`,
  `al` ("et al." citation), `multiphase` (the core domain term). None are
  misspellings. `inst/WORDLIST` lists these for the `spelling` unit test;
  the CRAN aspell check has no package-level suppression mechanism, so
  this NOTE recurs by design on every submission.

## Background

* Package contains five clinical reference datasets (~120 KB compressed)
  used in vignettes and parity tests against the original C/SAS HAZARD
  program.
* Vignettes use Quarto (`VignetteBuilder: quarto`).
* A subset of tests (SAS parity, multi-start multiphase fits, OMC
  raw-file derivation) are guarded with `skip_on_cran()` to stay within
  CRAN runtime budgets; the full suite runs on GitHub Actions.
* DOI URLs in vignettes (doi.org) may return HTTP 403 to automated
  checkers but resolve correctly in browsers.
