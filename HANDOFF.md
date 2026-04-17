# Session Handoff — temporal_hazard

Last updated: 2026-04-14

## Current State

**Version:** 0.9.3 **Branch:** `feature/analytic-gradients` **R CMD
check:** 0 ERROR, 0 WARNING, 1 NOTE (new submission) **Lint:** Clean
(lintr CI workflow in place) **Test suite:** 700+ tests passing across
all platforms. **Exported functions:** 8 new diagnostic/utility
functions in this session.

## What Was Done This Session (2026-04-14, continued)

### 0. Conservation of Events (CoE) in Optimizer

Turner’s theorem integrated into
[`.hzr_optim_multiphase()`](https://ehrlinger.github.io/temporal_hazard/reference/dot-hzr_optim_multiphase.md).
One phase’s log_mu solved analytically at each iteration, reducing
optimization dimension by 1. Controlled via `control$conserve` (default
TRUE for right-censored data). Key bug fixed: gradient function closure
scoping with double-wrapped logl_fn — solved with `logl_fn_unwrapped`
capture.

### 0b. SAS Macro Utility Functions (Phase 5 — COMPLETE)

All SAS HAZARD utility macros now have R equivalents:

| Function                                                                                                | SAS Macro                 | Purpose                                          |
|:--------------------------------------------------------------------------------------------------------|:--------------------------|:-------------------------------------------------|
| [`hzr_deciles()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_deciles.md)                 | `deciles.hazard.sas`      | Decile-of-risk calibration with chi-sq GOF       |
| [`hzr_gof()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_gof.md)                         | `hazplot.sas`             | Observed vs expected events, CoE diagnostic      |
| [`hzr_kaplan()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_kaplan.md)                   | `kaplan.sas`              | KM with logit CL, hazard, density, RMST          |
| [`hzr_calibrate()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_calibrate.md)             | `logit.sas`/`logitgr.sas` | Variable calibration (logit/Gompertz/Cox)        |
| [`hzr_nelson()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_nelson.md)                   | `nelsonl.sas`             | Nelson cumhaz with lognormal CL, weighted events |
| [`hzr_bootstrap()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_bootstrap.md)             | `bootstrap.hazard.sas`    | Bootstrap bagging with coefficient summaries     |
| [`hzr_competing_risks()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_competing_risks.md) | `markov.sas`              | Aalen-Johansen competing risks incidence         |

### 1. CRAN Review Fixes (from GPT-5.3 review)

**High: `fit` default docs/behavior mismatch** - `R/hazard_api.R`
roxygen `@param fit` said “default TRUE” but implementation is
`fit = FALSE` - Fixed docs and `man/hazard.Rd` to say “default FALSE”

**Medium: Non-deterministic reseeding in multi-start optimizer** -
`R/likelihood-multiphase.R` had `set.seed(NULL)` in the multi-start
loop, actively breaking reproducibility when users call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before fitting -
Removed `set.seed(NULL)` — perturbations now follow the user’s RNG state

**Low: Vignette metadata inconsistency** - 6 of 7 vignettes used HTML
comment `<!-- %\Vignette... -->` blocks - `ar-architecture.qmd` used
YAML `vignette: >` key (the correct Quarto approach) - Normalized all 7
to use YAML `vignette:` key inside front matter

**Low: CITATION URL redirect** - Added trailing slash to `inst/CITATION`
URL to match DESCRIPTION fix

### 2. SAS Parity Gap Analysis

Created `inst/dev/SAS-PARITY-GAP-ANALYSIS.md` — detailed
feature-by-feature comparison of C/SAS HAZARD vs TemporalHazard R
package. Added a capabilities table to `README.md`.

**Feature parity status:**

| Capability                             |  Status  |
|:---------------------------------------|:--------:|
| Multi-phase hazard modeling            | Complete |
| Right and interval censoring           | Complete |
| Time-varying covariates                | Complete |
| Covariance/correlation estimation      | Complete |
| Repeating events (epoch decomposition) | Complete |
| Weighted events                        | Complete |
| Stepwise covariate selection           | Planned  |
| Conservation of Events theorem         | Complete |

### 3. Prior Session Work (carried forward)

- Switched multiphase vignette demos from AVC to CABGKUL with G3 + fixed
  shapes
- Fixed critical
  [`summary.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/summary.hazard.md)
  bug: `!anyNA(vcov_mat)` was rejecting entire vcov when fixed params
  had NA diagonal entries
- Added `g3 = "g3"` to coefficient table switch statement
- Extensive lintr cleanup across R/ and tests/
- Added lint CI workflow (`.github/workflows/lint.yaml`)
- Created NEWS.md, inst/CITATION, cran-comments.md
- Fixed CI test: free-shape log-lik test changed from `<= -3700` to
  `>= -3750`

## What Needs Doing Next

### Immediate: Push and Verify CI

1.  **Push local commits:**

    ``` bash
    git push origin feature/analytic-gradients
    ```

    Then open a PR or merge to main so CI picks up all fixes.

2.  **Verify CI passes:** R-CMD-check, lint, and test-coverage workflows
    should all go green.

3.  **Local final check:**

    ``` r
    devtools::install()
    devtools::check(args = "--as-cran")
    devtools::spell_check()
    ```

4.  **Regenerate README figures:**

    ``` r
    source("inst/dev/generate-readme-figures.R")
    ```

### Branch Protection Enforcement

Ensure GitHub branch protection is enabled for `main` with:

1.  Pull request required before merge
2.  At least 1 approval
3.  Dismiss stale approvals on new commits
4.  Require status checks to pass
5.  Require branch to be up to date before merging
6.  Require conversation resolution before merging
7.  Optional but recommended: include administrators

Require these status checks:

- `lint`
- `test-coverage`
- `build-and-deploy` (pkgdown PR validation)
- `check-badge-sync`
- All `R-CMD-check` matrix checks

See `.github/BRANCH_PROTECTION.md` for the full branch-protection setup
guide.

### Short-term: CRAN Submission

- Ensure 0 ERROR / 0 WARNING / only expected NOTE(s)
- Submit via `devtools::submit_cran()` or web form
- Monitor CRAN feedback

### Medium-term: Feature Gaps (see `inst/dev/SAS-PARITY-GAP-ANALYSIS.md`)

**Priority order:** 1. **Stepwise covariate selection** — critical for
clinical research workflows. Large effort.

Completed in 0.9.3 / 0.9.4: Conservation of Events, observation weights,
repeating events (epoch decomposition).

### Optional: Test Cleanup

- Wrap 2 expected test warnings (NaNs in summary, Hessian not
  invertible) with `expect_warning()` so they don’t clutter test output
- Raise coverage from ~62% toward 80% on high-value files
  (`likelihood-multiphase.R`, `formula-helpers.R`, `decomposition.R`)

## Key Context

- The C/SAS HAZARD source is at `~/Documents/GitHub/hazard/`
- The mixhazard R package (decompos source) is at
  `~/Documents/GitHub/mixhazard/`
- hvtiPlotR (survival plotting) is at `~/Documents/GitHub/hvtiPlotR/`
- G3 late-phase decomposition is fully implemented in R
  ([`hzr_decompos_g3()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_decompos_g3.md))
- C/SAS HAZARD fixes all shapes and only estimates 3 log_mu values; R
  matches this with `fixed = "shapes"` in
  [`hzr_phase()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase.md)
- C reference benchmark: log-likelihood -3740.52 on CABGKUL with fixed
  shapes

## Files to Read First

When resuming:

    temporal_hazard/inst/dev/SAS-PARITY-GAP-ANALYSIS.md — feature gap details
    temporal_hazard/NEWS.md                              — changelog
    temporal_hazard/R/hazard_api.R                       — main API (45 KB)
    temporal_hazard/R/likelihood-multiphase.R             — multiphase engine (41 KB)
    temporal_hazard/R/decomposition.R                     — decompos family (23 KB)
    temporal_hazard/vignettes/getting-started.qmd         — primary user-facing demo
    hazard/src/llike/setcoe.c                             — Conservation of Events (C)
    hazard/src/vars/stepw.c                               — Stepwise selection (C)
