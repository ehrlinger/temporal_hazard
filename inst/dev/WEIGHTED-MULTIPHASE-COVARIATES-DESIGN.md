# Weighted Multiphase + Covariates — Coverage Design

**Date:** 2026-06-02
**Roadmap item:** Phase 7c — "Weighted multiphase + covariates" (DEVELOPMENT-PLAN.md)
**Type:** Test coverage (no new exported API)

## Problem

Every weighted multiphase test in `tests/testthat/test-weights.R` passes
`covariate_counts = c(early = 0L, constant = 0L)` and
`x_list = list(early = NULL, constant = NULL)` — i.e. **intercept-only**. No
weighted multiphase test exercises a covariate actually entering a phase, at
any level (likelihood, analytic gradient, or end-to-end fit). The analytic
gradient computes a per-phase `eta_j = log_mu + x_j %*% beta_j` and a
corresponding beta-score (`R/likelihood-multiphase.R:604`), so the weighted
beta-score path has never been verified.

This is the established pattern's last open multiphase coverage gap. Two of the
three coverage gaps closed in the v1.1.0 dev cycle harbored real bugs, so this
work may surface a latent defect rather than merely add a green test.

## Invariant

Integer-weight duplication parity, identical to the intercept-only tests:

    L(theta; D, w) == L(theta; expand(D, w), 1)

where `expand(D, w)` duplicates row `i` exactly `w_i` times
(`idx <- rep(seq_len(n), times = w)`). The design matrix is duplicated
row-for-row alongside `time`/`status`, so the invariant is exact for the
likelihood, gradient, MLE, and vcov.

Fractional / inverse-probability weights are **out of scope** — tracked
separately as roadmap item 7a.

## Fixtures

One multiphase setup reused across the new tests:

- Phases: `early = hzr_phase("cdf", ..., fixed = "shapes")` + `constant`.
- Covariates in **two** phases (one in `early`, one in `constant`) to also
  exercise the phase-specific-covariate path under weights.
- Integer weights `w in {1, 2, 3}`; `idx <- rep(seq_len(n), times = w)`
  duplicates `time`, `status`, and each phase's design matrix in lockstep.

## Tests

| # | Level | Check | File |
|---|---|---|---|
| 1 | Likelihood | `.hzr_logl_multiphase(weights = w, x_list, covariate_counts > 0)` == duplicated-row LL | `test-weights.R` |
| 2 | Gradient | `.hzr_gradient_multiphase` (weighted, covariates) == `numDeriv::grad` of the weighted LL | `test-weights.R` |
| 3 | Fit, CoE off | `hazard(Surv ~ x, weights = w, control = list(conserve = FALSE))` coef + vcov == duplicated-row fit | `test-weights.R` |
| 4 | Fit, CoE on | same fit with `conserve = TRUE` — weighted CoE + covariate interaction | `test-conservation-of-events.R` |

Tests 1–3 extend the multiphase section of `test-weights.R`. Test 4 lands in
`test-conservation-of-events.R`, where weighted-CoE parity already lives.

Tolerances and `skip_on_cran()` / `skip_if_not_installed("numDeriv")` guards
match the surrounding tests in each file.

## Out of scope

- Fractional / IPS weights (item 7a).
- `predict(decompose = TRUE, se.fit = TRUE)` weighted path (separate 7c row).
- Stepwise-under-weights and interval-censoring parity (separate 7c rows).

## Failure handling

If any test fails, it has surfaced a latent bug in the weighted-covariate
multiphase path. Stop and switch to systematic debugging; do not relax the
invariant or the tolerance to make a failing test pass.
