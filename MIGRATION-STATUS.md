# TemporalHazard Migration Status

Last updated: 2026-04-07

This file is the single source of truth for migration progress in this repository.

## Executive Summary

The pure-R migration is through M4 with core features, testing, and documentation substantially complete.

- M0 complete: package scaffold, CI, documentation tooling, math primitives
- M1 complete: public constructor/predict/print API, argument mapping, migration vignette baseline
- M2 complete: Weibull fitting core, analytical gradient, optimizer path, prediction expansion, synthetic golden fixtures, parity harness
- M3 complete: exponential, log-logistic, and log-normal support; mixed censoring; piecewise time-varying coefficients; expanded parity and edge-case coverage; legacy binary parity helpers
- M4 complete: `summary.hazard()` method, formula interface (`Surv(...) ~ predictors`), comprehensive vignettes, release documentation, platform-aware CI

Current validation state:

- 275 tests passing; 1 platform-aware skip (Unix legacy binary parity)
- `R CMD check --no-manual` status: OK
- Working tree: clean at the time this file was updated

## What Is Done

### Core modeling

- Weibull likelihood, gradient, and optimizer are implemented
- Exponential likelihood, gradient, and optimizer are implemented
- Log-logistic likelihood, gradient, and optimizer are implemented
- Log-normal likelihood, gradient, and optimizer are implemented
- `hazard()` supports fitting across the implemented distributions

### Prediction support

- `predict.hazard()` supports `linear_predictor`
- `predict.hazard()` supports `hazard`
- `predict.hazard()` supports `survival`
- `predict.hazard()` supports `cumulative_hazard`
- Piecewise time-varying coefficient support is implemented in fit and predict paths

### Censoring support

- Right censoring supported
- Left censoring supported
- Interval censoring supported
- Mixed censoring scenarios covered across Weibull, exponential, log-logistic, and log-normal

### Parity and fixtures

- Synthetic golden fixtures are present in `inst/fixtures/`
- Core parity tests are present for refit reproducibility and prediction consistency
- Legacy binary helpers exist in `R/parity-helpers.R`
- Legacy execution now matches shell-redirection usage:
  - `hazard.exe < prefix.sas > prefix.lst 2>&1`
  - `hazpred.exe < prefix.sas > prefix.lst 2>&1`
- Generated `%HAZARD(...)` and `%HAZPRED(...)` control files are implemented for parity workflows
- Legacy binary tests skip cleanly when the executable is unavailable or non-runnable

### Public API already present

- `hazard()`
- `predict.hazard()`
- `print.hazard()`
- `summary.hazard()`
- `coef.hazard()`
- `vcov.hazard()`
- `hzr_argument_mapping()`
- Formula interface: `hazard(Surv(time, status) ~ predictors, data = df, ...)`

## What Is Not Done Yet

These are the remaining migration items that still look materially open.

### M4 API stabilization

- Formula interface is implemented: `hazard(Surv(time, status) ~ predictors, data = df, ...)`
- Public API freeze has not been documented yet
- `confint.hazard()` is still optional and not implemented

### Documentation and release polish

- The repository had conflicting milestone status documents; this file resolves that going forward
- Vignettes exist as Quarto files, but they still need a final pass to confirm the migration narrative and examples are complete and consistent
- Full pkgdown reference and worked-example review remains open
- `R CMD check --as-cran` is still a separate final-release gate and is not recorded here as complete

### Ongoing parity caution

- The parity infrastructure is in place and tested, but broader validation against real legacy executables and production-style datasets is still worth doing where those binaries are available

## Milestone Reconciliation

### M0

Complete.

Included scaffold, CI, pkgdown setup, roxygen workflow, and numerical helper primitives.

### M1

Effectively complete for the original migration-baseline scope.

The only items that were historically tracked here and still matter are operational/project-hosting tasks, not package-core migration work.

### M2

Complete.

The older M2 note is obsolete because it still describes gradients, fixtures, and parity as pending even though they are implemented.

### M3

Complete.

The older M3 note is directionally correct on scope, but it understates the current validation count and predates the latest legacy parity work.

### M4

Current active milestone.

Recommended focus:

1. Add formula interface support
2. Review and tighten vignette/example completeness
3. Freeze and document the supported public API
4. Run final release-oriented checks, including `--as-cran`

## Source Files For Current State

- `R/hazard_api.R`
- `R/likelihood-weibull.R`
- `R/likelihood-exponential.R`
- `R/likelihood-loglogistic.R`
- `R/likelihood-lognormal.R`
- `R/parity-helpers.R`
- `tests/testthat/`
- `inst/fixtures/`

## Archived Documents

The following files are retained as historical notes only:

- `TODO.md`
- `M2-STATUS.md`
- `M3-STATUS.md`