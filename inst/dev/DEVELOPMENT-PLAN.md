# TemporalHazard Development Plan

Last updated: 2026-05-13

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

### 4a. Conservation of Events Theorem — COMPLETE (v0.9.3)

Turner's theorem integrated into `.hzr_optim_multiphase()`. One phase's
log_mu solved analytically at each iteration via `.hzr_conserve_events()`.
Controlled via `control$conserve` (default TRUE). Currently supports
right-censored + exact-event data; interval/left-censored extension is a
follow-up item. Key functions: `.hzr_conserve_events()`,
`.hzr_select_fixmu_phase()`, `.hzr_log_mu_positions()`.

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

### 4c. Observation Weights -- COMPLETE (v0.9.6)

`weights` argument on `hazard()` applies Fisher weighting to the
log-likelihood across all five distributions (Weibull, exponential,
log-logistic, log-normal, multiphase). Threaded through every
distribution's likelihood and analytic gradient, the multiphase
Conservation-of-Events adjustment (`.hzr_conserve_events()` weights
per-phase cumhaz sums so Turner's adjustment is on the same scale as
the weighted event count), and the numerical-gradient fallbacks for
mixed censoring. Integer weights reproduce the row-duplicated fit to
optimizer tolerance for every distribution.

### 4d. Repeating Events (Epoch Decomposition) -- COMPLETE (v0.9.7)

`Surv(start, stop, event)` is parsed (start -> time_lower, stop -> time,
event -> status) and scored as `H(stop) - H(start)` per epoch in the
Weibull and multiphase likelihoods. The trivial `Surv(0, t, d)` case
degenerates to `H(stop)` and recovers the plain-Surv fit exactly. Full
wire-up is tracked under 4f below.

### 4e. Complete `weights` wire-up (exp / log-logistic / log-normal + CoE) -- COMPLETE (v0.9.6)

Every LL term in `R/likelihood-exponential.R`,
`R/likelihood-loglogistic.R`, and `R/likelihood-lognormal.R` now
multiplies by `weights[idx_*]`; analytic gradients use weighted
building blocks (`w * status`, `w * cumhaz`, `w * (1+delta) * pw`
for log-logistic, `w * (1-delta) * mills` for log-normal). Numerical
gradient fallbacks forward `weights` through. The
`dist %in% c("weibull", "multiphase")` guard in `hazard()` has been
removed.

`.hzr_conserve_events()` and `.hzr_select_fixmu_phase()` take an
optional `weights` argument applied to per-phase cumhaz sums; the
`all(weights == 1)` guard on CoE in `R/likelihood-multiphase.R` is
gone. Duplication-parity tests for each distribution are in
`test-weights.R`; weighted-CoE parity (CoE on and off) lives in
`test-conservation-of-events.R`.

### 4f. Complete repeating-events wire-up (counting-process LL) -- COMPLETE (v0.9.7)

`.hzr_logl_weibull()` and `.hzr_logl_multiphase()` apply
`H(stop) - H(start)` to event and right-censored terms via a
`cumhaz_start` term that is zero when no `time_lower` is supplied
(plain-Surv path). Analytic gradients for both distributions add the
matching `-d H(start)/d theta` contributions (Weibull closed-form;
multiphase per-phase `Phi_j(start)` and shape derivatives, reusing
the existing G3 finite-difference machinery). The `hazard()` guard on
`start > 0` is removed; the three split-invariance tests in
`test-repeating-events.R` are live.

Remaining follow-ups (low priority): vignette section on
epoch-decomposed longitudinal data, and a SAS-parity fixture for a
counting-process reference fit.

### 4g. Prediction confidence limits (delta method) -- COMPLETE (v0.9.8)

`predict.hazard()` takes `se.fit = FALSE, level = 0.95`.  When
`se.fit = TRUE`, the return value is a data frame with columns
`fit`, `se.fit`, `lower`, `upper`.  Per-row SEs come from the delta
method: `se_i = sqrt( (J V J^T)_{ii} )` where `J[i, ] = d eta_i / d theta`
and `V = vcov(theta_hat)`.

* **Jacobians.** Closed-form for Weibull
  (`.hzr_predict_jacobian_weibull()`) and multiphase
  (`.hzr_predict_jacobian_multiphase()`, reusing
  `.hzr_phase_derivatives()` / `.hzr_g3_phase_derivatives()`);
  `numDeriv::jacobian()` fallback for exponential / log-logistic /
  log-normal via a per-call cumhaz closure.
* **Scale.** `hazard` and `cumulative_hazard` use log-scale CLs;
  `survival` uses `log(-log(S)) = log H` CLs so 0 <= lower <= upper
  <= 1; `linear_predictor` is symmetric on the natural scale.  Matches
  SAS `hzp_calc_haz_CL.c` / `hzp_calc_srv_CL.c`.
* **Fixed parameters.** `fixed = "shapes"` and CoE leave NA rows/cols
  in the expanded vcov.  The sandwich restricts to the free-parameter
  submatrix and zero-indexes the corresponding Jacobian columns, so
  CLs reflect only free-parameter uncertainty.
* **Backward compat.** `se.fit = FALSE` preserves the pre-0.9.8
  scalar-vector return shape.  `decompose = TRUE` with `se.fit = TRUE`
  is rejected with a clean error.

### 4b scheduling note

With 4a, 4b, 4c, 4d, 4e, 4f, and 4g complete, the remaining SAS-parity
gaps are: density / quantile prediction types, `OUTEST` / `OUTVCOV`
dataset export, `anova.hazard()` S3 method, `NOCOV` / `NOCOR` summary
flags, and the stepwise enhancements (FAST screening, multi-step MOVE
trace) tracked separately.

---

## Phase 5: Utility Functions (SAS Macro Equivalents) — COMPLETE (v0.9.3)

The SAS HAZARD system ships utility macros for non-parametric estimation,
goodness-of-fit, variable calibration, bootstrap inference, and competing
risks. These cover the full analytical workflow around the core fitting
engine. R equivalents should be exported functions in the `hzr_` namespace.

### SAS macro inventory and R mapping

| SAS Macro | Purpose | R Equivalent | Priority |
|:---|:---|:---|:---:|
| `kaplan.sas` | KM survival with exact logit CL | `hzr_kaplan()` | Done |
| `hazplot.sas` | Observed vs predicted overlay + CoE GOF | `hzr_gof()` | Done |
| `deciles.hazard.sas` | Calibration by predicted survival deciles | `hzr_deciles()` | Done |
| `chisqgf.sas` | Chi-square GOF test (obs vs expected) | (part of `hzr_deciles()`) | Done |
| `nelsonl.sas` | Nelson estimator for weighted/repeated events | `hzr_nelson()` | Done |
| `nelsont.sas` | Nelson estimator for terminating events | `hzr_nelson()` (unified) | Done |
| `bootstrap.hazard.sas` | Bootstrap bagging with variable selection | `hzr_bootstrap()` | Done |
| `bootstrap.summary.sas` | Summarize bootstrap coefficient distributions | (part of `hzr_bootstrap()`) | Done |
| `logit.sas` | Calibrate continuous variables (decile logit) | `hzr_calibrate()` | Done |
| `logitgr.sas` | Calibrate by natural grouping variable | `hzr_calibrate()` (with `by=`) | Done |
| `markov.sas` | Competing risks cumulative incidence (Greenwood) | `hzr_competing_risks()` | Done |
| `plot.sas` | Publication-quality SAS plots | Not needed (ggplot2) | — |

### 5a. `hzr_kaplan()` — Kaplan-Meier with Exact Confidence Limits

**SAS reference:** `kaplan.sas`

R's `survival::survfit()` provides KM estimates, but the SAS macro adds
logit-transformed exact confidence limits (more accurate than Greenwood
in the tails), hazard rate estimation at each interval, life integral
computation, and stratified log-rank testing. A thin wrapper around
`survfit()` that adds these extras and returns a tidy data frame would
match the SAS output structure used by `hazplot.sas` and `deciles.hazard.sas`.

**Effort:** Small — mostly formatting `survfit()` output with added CL
transform.

### 5b. `hzr_gof()` — Goodness-of-Fit: Observed vs Predicted

**SAS reference:** `hazplot.sas`

The central validation tool in the HAZARD workflow. Given a fitted `hazard`
object, compute at each observed event time:
- Parametric predicted survival, cumulative hazard, and hazard rate
- Nonparametric (KM or Nelson) estimates at the same times
- Cumulative observed vs expected event counts (conservation of events)
- Chi-square GOF p-value for the obs/exp comparison
- Per-phase hazard component percentages

Returns a data frame suitable for plotting parametric vs nonparametric
overlays with `ggplot2`. Optionally supports stratified GOF by a grouping
variable.

**Effort:** Medium — requires combining `predict.hazard()` output with
`hzr_kaplan()` output and implementing the CoE event-counting algorithm.

### 5c. `hzr_deciles()` — Decile-of-Risk Calibration

**SAS reference:** `deciles.hazard.sas`

Given a fitted model and a time point, rank patients into deciles by
predicted survival, then compare actual vs expected event counts per
decile with chi-square GOF. The inference vignette already demonstrates
this workflow manually — wrapping it as a function standardizes the
output and enables programmatic calibration checks.

**Effort:** Small — the logic is already demonstrated in the vignette;
needs packaging as an exported function.

### 5d. `hzr_nelson()` — Nelson Cumulative Hazard Estimator

**SAS reference:** `nelsonl.sas`, `nelsont.sas`

Wayne Nelson's estimator for cumulative hazard with lognormal confidence
limits. Handles both terminating (single) events and weighted/repeated
events. Not available in base `survival` — `survfit()` uses the
Breslow estimator which differs in variance estimation. Important for
validating fits on weighted or repeating-event data.

**Effort:** Medium — the Nelson algorithm is straightforward but the
lognormal CI transform and weighted-event handling need careful
implementation.

### 5e. `hzr_bootstrap()` — Bootstrap Inference

**SAS reference:** `bootstrap.hazard.sas`, `bootstrap.summary.sas`

Resample data with replacement, refit the hazard model on each replicate,
and accumulate coefficient distributions. Supports fractional sampling
(bagging), time-varying covariate data (resample by patient ID), and
optional stepwise selection per replicate. Returns a tidy data frame of
per-replicate coefficients with summary statistics (frequency of selection,
mean, SD, CI).

**Effort:** Medium — the resampling loop is simple, but handling TVC data,
failed fits, and optional stepwise selection adds complexity. The inference
vignette already shows manual bootstrap; this wraps it as a function.

### 5f. `hzr_calibrate()` — Variable Calibration

**SAS reference:** `logit.sas`, `logitgr.sas`

Group a continuous covariate into quantile bins, compute event probability
per bin, and apply logit (or Gompertz/Cox) transform. Used for assessing
functional form before model entry — does the covariate-outcome
relationship look linear on the logit scale, or does it need a transform?
Supports stratification by a grouping variable.

**Effort:** Small — simple rank + aggregate + transform pipeline.

### 5g. `hzr_competing_risks()` — Competing Risks Incidence

**SAS reference:** `markov.sas`

Compute cumulative incidence of competing events using Greenwood's variance
formula. Unlike 1-KM (which overestimates incidence when competing risks
exist), this provides correct marginal incidence for each event type with
standard errors. The valves dataset (death, PVE, reoperation) is a natural
test case.

**Effort:** Medium — the Greenwood variance matrix computation is
moderately complex. Could also wrap `cmprsk::cuminc()` if a dependency
is acceptable.

---

## Phase 6: Documentation Gaps — PLANNED

The SAS HAZARD documentation walks users through a disciplined analytical
sequence. The R vignettes cover the same pieces but have gaps relative to
what a SAS HAZARD veteran would expect.

### 6a. Complete Clinical Analysis Walkthrough

**Priority:** High

A new vignette that mirrors the SAS `ac.*` → `hz.*` → `hm.*` → `hp.*` →
`hs.*` workflow on one dataset (AVC):

1. Nonparametric baseline (KM life table)
2. Shape fitting: start simple (Weibull), build to multiphase, show how
   to fix overparameterized shapes
3. Variable selection: manual screening with `hzr_calibrate()` and
   univariable logistic, then multivariable model building
4. Prediction: patient-specific risk profiles, covariate sensitivity
5. Validation: `hzr_gof()` overlay, `hzr_deciles()` calibration,
   conservation-of-events check

This is the single most impactful documentation addition — it teaches
the analytical discipline, not just the API.

### 6b. Cox Regression Comparison

**Priority:** Medium

A section (in getting-started or a standalone vignette) making the case
for multiphase parametric modeling vs Cox regression:

| Feature | Cox (semi-parametric) | TemporalHazard (parametric) |
|---|---|---|
| Proportional hazards | Required | Not required |
| Baseline hazard | Unspecified | Fully parametric |
| Multiple hazard phases | Not supported | Supported |
| Patient-specific prediction | Approximate | Exact |
| Event conservation check | Not available | Planned |
| Interval censoring | Not standard | Supported |

The SAS introduction makes this comparison prominently. R users coming
from `coxph()` need to understand why they'd switch.

### 6c. Convergence Troubleshooting

**Priority:** Medium

A section in the fitting vignette or math foundations covering:
- How to choose starting values from the KM cumulative hazard plot
- When to fix shape parameters (`fixed = "shapes"`) vs estimate freely
- Multi-start strategy and `control$n_starts`
- Signs of overparameterization (boundary estimates, huge SEs, NaN)
- Optimization method differences (the SAS docs walk through Newton vs
  Quasi-Newton vs Steepest Descent)

### 6d. Interval Censoring and Left-Truncation Examples

**Priority:** Low

Worked examples demonstrating interval-censored and left-truncated data
with `Surv(time1, time2, status)` syntax. The math foundations vignette
covers the theory but there are no hands-on examples.

---

## Phase 7: Production Model Parity & Edge Cases — PLANNED

Post-CRAN hardening against real production models from the original C/SAS
HAZARD system. Goal: high confidence in both codebases for the model types
that appear in clinical production, not just the curated examples shipped
with the C source.

> **Cross-codebase note:** These items are logged in both this document and
> `docs/PRODUCTION-TEST-PLAN.md` in the `hazard` C/SAS repo. A single `.lst`
> reference output validates both codebases simultaneously. CCF-local fixtures
> live at `tests/production/` (gitignored) in the hazard repo.

### 7a. Rajeswaran / Blackstone Production Model Parity

**Contact:** Rajeswaran (Rajes) at CCF — has the "hard" production hazard
models run against the C HAZARD binary in clinical analyses.

**Goal:** Collect 3–5 production `.lst` outputs from Blackstone's
Rajeswaran-authored models; port each to `hazard()` calls; add
`skip_on_cran()` parity tests with tolerances matching the existing suite.

Candidate model types not yet exercised:
- 4-phase models (early + early2 + constant + late)
- Phase-specific covariate models with multiple covariate phases
- Models with shaping modifiers (`/S`, `/I`, `/E`) — currently deferred from
  primary parity (see PRE-CRAN-PARITY-INVENTORY.md gap §1)
- High-dimensional covariate fits (12+ free parameters)

**Action:** Schedule meeting with Rajes; request SAS driver scripts + `.lst`
reference outputs for 3–5 candidate models from production library.

**Sequencing dependency:** Collect Rajes examples *before* installing
`hazard` v4.4.6 at CCF. Once the CCF production install is live, run the
new driver scripts against the installed binary to generate fresh `.lst`
captures; add those as `skip_on_cran()` parity fixtures. The CCF install
is blocked on this collection step.

### 7b. Weighted Events — Additional Parity Coverage

The OMC dataset (`hz.tm123.OMC`) exercises case weights via the PRIMISOL
morbidity endpoint, and `4c`/`4e` verified integer-weight duplication parity.
Additional coverage needed to close confidence gaps:

- **Fractional weights** (non-integer, e.g. inverse-probability weights): add
  a parity test against a SAS run with `WEIGHT` statement using fractional values
- **Weighted multiphase with covariates**: current weighted tests are intercept-only;
  need a covariate fit with weights to confirm the weighted score/gradient is
  correct under the full CoE + covariate path
- **Weighted competing risks** (`hzr_competing_risks()` with case weights):
  not currently tested with non-unit weights
- **Source:** Look for a CCF production model that uses fractional IPS weights
  or aggregated-cell weights (e.g. from a registry with pre-computed weights
  rather than individual patient records)

### 7c. New Example Datasets

**Hsich / UNOS (post-heart-transplant mortality)**
- 3-phase early + constant + late structure
- SRTR Standard Analysis Files (requires DUA — ship cleaned summary or
  synthetic derivative in `data/`, keep raw pull off-repo)
- Eileen Hsich (CCF) is the contact; check for published paper with
  parameter estimates to anchor reproducibility

**Sargent CABG mortality paper**
- Published analysis of death after CABG using multiphase hazard model
- If SAS HAZARD was used, published MLEs serve as reference values
- Actionable: locate paper, extract published LL + MLE table, build
  vignette that reproduces the table with `hazard()`
- Distinct from KUL (Belgium/Leuven) which is already the canonical
  3-phase example; Sargent adds a second independent CABG validation

**Acute aortic dissection — 4-phase model**
- First package example exercising a 4-phase decomposition end-to-end
- Clinical phase structure:
  1. **Operative** (G1/CDF): very high immediate hazard, hours to days
  2. **Early** (G1/CDF): subacute complications, weeks to ~90 days
  3. **Constant** (G2): background mortality plateau
  4. **Late** (G3): structural deterioration, false-lumen expansion,
     reoperation risk rising with time
- Data source candidates: IRAD (International Registry of Acute Aortic
  Dissection) published series; large single-center Type A repair
  publications with parameter estimates; or CCF internal series
- Value: first real test of 4-phase optimizer convergence, CoE algebra
  with 4 phases, and `predict(..., decompose = TRUE)` output at scale
- Actionable: identify published paper with LL + MLE table, or request
  SAS driver from Blackstone/CCF production library

### 7d. Code Areas Requiring Deeper Examination

The following components are implemented but have limited parity coverage
or known numerical gaps. Priority order for investigation:

| Area | Risk | Status |
|---|---|---|
| **4-phase CoE algebra** | `.hzr_select_fixmu_phase()` selects one constrained phase from N; behavior with N=4 untested | No test |
| **`hzr_decompos()` near-boundary cases** | Sign dispatch has 6 cases; edge values (nu≈0, m≈0, t_half→0, t_half→∞) may fall between branches | Partial — Case 2L verified |
| **Phase-specific covariates** | `hzr_phase(formula = ~)` ships but no SAS parity fixture exercises it | No parity test |
| **`hzr_stepwise()` on multiphase** | Phase-specific variable entry path — no parity test against SAS SELECTION | No parity test |
| **`hzr_bootstrap()` with weights** | Non-unit weights path through bootstrap resampling untested | No test |
| **`predict(..., decompose=TRUE, se.fit=TRUE)`** | Currently blocked with clean error — delta-method Jacobian needs per-phase extension | Blocked |
| **Interval censoring under multiphase** | Code path exists; no real-data or SAS parity fixture | No parity test |
| **Hessian stability at 12+ parameters** | Numerical Hessian inversion can become ill-conditioned; `hm.death.AVC.deciles` (13 params) passes but borderline | Passing, fragile |
| **`hzr_competing_risks()` with weights** | Greenwood variance with case weights not tested | No test |
| **Weighted multiphase + covariates** | All weighted parity tests are intercept-only; covariate + weight combination untested | No test |

---

## Phase 8: Performance and Extensions — FUTURE

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

## Cross-codebase documents

| Document | Location | Purpose |
|:---|:---|:---|
| Production test plan | `hazard/docs/PRODUCTION-TEST-PLAN.md` | CCF-local fixture collection plan; shared `.lst` validates both codebases |
| hazard corpus acceptance | `hazard/docs/VALIDATION_PLAN.md` | Public example corpus + black-box harness |
| hazard v4.4.6 plan | `hazard/docs/superpowers/plans/2026-04-29-v4.4.6-implementation-plan.md` | v4.4.6 implementation tasks |
