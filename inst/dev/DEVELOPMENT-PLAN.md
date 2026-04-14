# TemporalHazard Development Plan

Last updated: 2026-04-14

This is the single source of truth for the development roadmap of the
TemporalHazard R package — a pure-R implementation of the multiphase
parametric hazard model of Blackstone, Naftel, and Turner (1986).

For detailed feature-level parity analysis against the original C/SAS HAZARD
system, see `SAS-PARITY-GAP-ANALYSIS.md` in this directory.

---

## Phase 1: Migration from C/SAS — COMPLETE

Ported the core modeling engine from compiled C code and SAS macro wrappers
to a pure-R package with modern idioms (S3 classes, formula interface,
roxygen2, testthat).

### Milestones

| Milestone | Scope | Status |
|:---|:---|:---:|
| M0 | Package scaffold, CI, pkgdown, math primitives | Complete |
| M1 | Public constructor/predict/print API, argument mapping, migration vignette | Complete |
| M2 | Weibull fitting, analytical gradient, optimizer, golden fixtures, parity harness | Complete |
| M3 | Exponential, log-logistic, log-normal; mixed censoring; time-varying covariates; parity helpers | Complete |
| M4 | `summary.hazard()`, formula interface (`Surv() ~ predictors`), vignettes, release docs, CI | Complete |

### Delivered capabilities

- 5 distributions: Weibull, exponential, log-logistic, log-normal, multiphase
- Right, left, and interval censoring with correct likelihood contributions
- Piecewise time-varying coefficients via `time_windows`
- Formula interface: `hazard(Surv(time, status) ~ age + nyha, data = df, ...)`
- S3 methods: `coef()`, `vcov()`, `summary()`, `predict()`, `print()`
- `hzr_argument_mapping()` for SAS/C parameter name translation
- Golden fixture regression tests and C binary parity infrastructure
- 7 Quarto vignettes covering getting started through mathematical foundations

**Detail:** `MIGRATION-STATUS.md`

---

## Phase 2: Multiphase Hazard Models — COMPLETE

Implemented the N-phase additive hazard decomposition that is the
distinguishing feature of the Blackstone/Naftel/Turner framework.

### Architecture

Additive cumulative hazard across J phases:

    H(t|x) = Sigma_j  mu_j(x) * Phi_j(t; t_half_j, nu_j, m_j)

where mu_j(x) = exp(alpha_j + x * beta_j) and Phi_j is the temporal shape
function for phase j.

### Delivered capabilities

- Generalized temporal decomposition family (`hzr_decompos()`) — all 6
  mathematical cases based on signs of nu and m
- Phase types: G1/CDF (early), G2/constant, G3 (late, unbounded)
- `hzr_phase()` constructor with `fixed = "shapes"` matching SAS `FIXparm`
- Phase-specific covariate formulas
- Multi-start optimizer with Nelder-Mead warm-up for fixed-shape problems
- `predict(..., decompose = TRUE)` for per-phase cumulative hazard
- Phase-grouped coefficient table in `summary.hazard()`
- 540+ tests, C binary parity validation on CABGKUL dataset
  (log-likelihood -3740.52 matching C reference with fixed shapes)

### Files

| File | Description |
|:---|:---|
| `R/decomposition.R` | `hzr_decompos()`, `hzr_decompos_g3()`, phase cumhaz/hazard |
| `R/phase-spec.R` | `hzr_phase()` constructor and helpers |
| `R/likelihood-multiphase.R` | Likelihood, gradient, multi-start optimizer |
| `R/hazard_api.R` | Multiphase dispatch, predict decompose, summary |

**Detail:** `MULTIPHASE_PLAN.md`

---

## Phase 3: CRAN Release — IN PROGRESS

Preparing the package for initial CRAN submission.

### Completed

- `R CMD check --as-cran`: 0 ERROR, 0 WARNING, 1 NOTE (new submission)
- lintr clean with CI enforcement (`.github/workflows/lint.yaml`)
- NEWS.md, inst/CITATION, cran-comments.md
- Fixed `fit` default docs/behavior mismatch
- Removed non-deterministic `set.seed(NULL)` in multi-start optimizer
- Normalized vignette metadata (YAML `vignette:` key across all 7 files)
- Fixed CITATION and DESCRIPTION URL redirects (trailing slash)

### Remaining

- Push local commits and verify all CI workflows pass green
- `devtools::spell_check()` — fix any typos
- Regenerate README figures: `source("inst/dev/generate-readme-figures.R")`
- Submit via `devtools::submit_cran()` or web form
- Monitor and respond to CRAN feedback

---

## Phase 4: SAS Feature Parity — PLANNED

Close the remaining feature gaps between the C/SAS HAZARD system and the
R package. See `SAS-PARITY-GAP-ANALYSIS.md` for detailed implementation
notes on each item.

### 4a. Conservation of Events Theorem

**Priority:** High
**Effort:** Medium-large
**Impact:** Improved numerical stability and convergence for all users

The CoE theorem (Turner) provides an explicit closed-form solution for one
scaling parameter (mu) given the other parameters, reducing the optimization
dimension by one. The C implementation lives in `setcoe.c` and `consrv.c`.

Tasks:
1. Implement CoE constraint solver: given shape params and data, solve for
   log_mu satisfying Sum(events) = Sum(CF(t_i))
2. Integrate into multiphase optimizer loop (solve mu analytically at each
   iteration before evaluating likelihood on reduced parameter set)
3. Extend to interval-censored and weighted data
4. Add `conserve` control option (default TRUE) matching SAS
   `CONSERVE/NOCONSERVE`
5. Validate against SAS output

### 4b. Stepwise Covariate Selection

**Priority:** High
**Effort:** Large
**Impact:** Critical for exploratory clinical research workflows

The C/SAS system has forward, backward, and two-way stepwise selection with
configurable SLENTRY/SLSTAY thresholds, MOVE limits, forced entry ordering,
and FAST approximate Wald screening. Selection operates within the hazard
model, accounting for multiphase structure.

Tasks:
1. Wald test infrastructure for individual covariate effects within
   multiphase framework
2. Forward step: fit with each candidate added, select smallest p-value
   below SLENTRY
3. Backward step: test each included variable, remove largest p-value above
   SLSTAY
4. Stepping loop with MOVE limits and convergence detection
5. Phase-specific variable selection (variables enter different phases
   independently)
6. FAST screening mode using approximate Wald updates (optimization)
7. User-facing API: `hzr_stepwise()` or `step.hazard()` method
8. Selection trace printing/reporting

### 4c. Observation Weights

**Priority:** Medium
**Effort:** Medium
**Impact:** Enables severity-weighted analyses

Tasks:
1. Add `weights` argument to `hazard()` (numeric vector)
2. Thread through all distribution-specific log-likelihood functions
3. Thread through gradient functions
4. Update Hessian computation
5. Integrate with CoE theorem extension
6. Validate against SAS output

### 4d. Repeating Events (Epoch Decomposition)

**Priority:** Medium
**Effort:** Small-medium
**Impact:** Enables longitudinal recurrent event analyses

The likelihood math for start-stop intervals is already implemented via
interval censoring. The main work is data interface and documentation.

Tasks:
1. Ensure `Surv(stime, time, status)` start-stop notation is correctly
   parsed and routed to interval-censoring likelihood path
2. Document the epoch-decomposition data preparation workflow
3. Add a vignette or example demonstrating repeating events
4. Validate against SAS output on a repeating-events dataset

---

## Phase 5: Performance and Extensions — FUTURE

Items that would improve the package beyond SAS parity.

- **Rcpp acceleration** for likelihood hot paths (decomposition, cumulative
  hazard evaluation) — currently pure R
- **Analytical gradients for multiphase** — currently uses numerical
  gradients via `numDeriv`; analytical score vectors would improve optimizer
  precision and speed
- **`confint.hazard()`** — profile likelihood or Wald-based confidence
  intervals as an S3 method
- **Residuals** — deviance, martingale, or Schoenfeld residuals for
  model diagnostics
- **Frailty / random effects** — shared frailty models for clustered data
  (distinct from the epoch-decomposition approach to repeating events)
- **Parallel bootstrap** — leverage `future`/`furrr` for bootstrap CI
  computation

---

## Related Documents

| File | Purpose |
|:---|:---|
| `inst/dev/MIGRATION-STATUS.md` | Phase 1 milestone detail and reconciliation |
| `inst/dev/MULTIPHASE_PLAN.md` | Phase 2 architecture, math, API design, build order |
| `inst/dev/SAS-PARITY-GAP-ANALYSIS.md` | Phase 4 feature-by-feature gap analysis with C/SAS code references |
| `HANDOFF.md` | Session-to-session context transfer (top-level, .Rbuildignored) |
| `NEWS.md` | User-facing changelog |
| `cran-comments.md` | CRAN submission notes |

## Reference Repositories

| Repo | Path | Purpose |
|:---|:---|:---|
| hazard | `~/Documents/GitHub/hazard/` | Original C/SAS source code |
| mixhazard | `~/Documents/GitHub/mixhazard/` | Generalized decompos() R source |
| hvtiPlotR | `~/Documents/GitHub/hvtiPlotR/` | Survival plotting companion package |
