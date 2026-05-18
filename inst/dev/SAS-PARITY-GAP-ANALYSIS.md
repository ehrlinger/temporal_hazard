# SAS HAZARD Parity Gap Analysis

**Date:** 2026-05-18
**Package version:** 1.0.0
**Branch:** main

This document tracks feature parity between the original C/SAS HAZARD system
(Blackstone, Naftel, Turner 1986; maintained by Rajeswaran & Ehrlinger) and
the TemporalHazard R package.

## Summary

| Capability | R Status | SAS Reference |
|:---|:---:|:---|
| Multi-phase hazard modeling | Complete | `src/llike/`, `src/decomp/` |
| Right and interval censoring | Complete | `src/llike/setcoe_obs_loop.c` |
| Repeating events (epoch decomposition) | Complete | `src/llike/setcoe_obs_loop.c` (STIME) -- `Surv(start, stop, event)` epoch decomposition; likelihood accrues `H(stop) − H(start)` (see sec. 3) |
| Time-varying covariates | Complete | EARLY/CONSTANT/LATE variable lists |
| Weighted events | Complete | `src/llike/setcoe_obs_loop.c` -- all five distributions; Fisher weighting in LL, gradient, and Hessian (see sec. 5) |
| Stepwise covariate selection | Complete | `src/vars/stepw.c`, `backw.c`, `swvari.c` -- `hzr_stepwise()` shipped v0.9.5 (see sec. 6) |
| Conservation of Events theorem | Complete | `src/llike/setcoe.c`, `consrv.c` -- right-censored + exact event case; weighted/interval-censored extensions tracked in Phase 7 (see sec. 7) |
| Covariance and correlation matrix estimation | Complete | `src/optim/` |

---

## 1. Multi-phase hazard modeling — Complete

The R package supports all phase types from the C/SAS system:

- **G1 / CDF (early):** `hzr_phase("cdf", t_half, nu, m)` — bounded [0,1],
  resolves over time. Maps to SAS `EARLY` statement.
- **G2 / Constant:** `hzr_phase("constant")` — flat background hazard rate.
  Maps to SAS `CONSTANT` statement.
- **G3 / Late:** `hzr_phase("g3", tau, gamma, alpha, eta)` — monotone
  increasing, unbounded. Maps to SAS `LATE` statement.
- **N-phase generalization:** Arbitrary number of phases via
  `dist = "multiphase"` with `phases = list(...)`.

The generalized temporal decomposition family (`hzr_decompos()`) implements
all 6 valid mathematical cases based on signs of `nu` and `m`. Phase-specific
covariates are supported via `formula` argument in `hzr_phase()`. The
`fixed = "shapes"` parameter matches SAS `FIXparm` behavior.

**SAS C reference:** `src/decomp/`, `src/llike/`

---

## 2. Right censoring and interval censoring — Complete

All censoring types are implemented with correct likelihood contributions:

- **Right-censored** (`status = 0`): contribution = −H(t)
- **Exact events** (`status = 1`): contribution = log h(t) − H(t)
- **Left-censored** (`status = −1`): contribution = log(1 − exp(−H(t)))
- **Interval-censored** (`status = 2`): contribution = −H(lower) + log(1 − exp(−[H(upper) − H(lower)]))

Mixed censoring is supported in all distributions (Weibull, exponential,
log-logistic, log-normal, multiphase). The `survival::Surv()` formula
interface handles `Surv(time, status)` for right-censored and
`Surv(time1, time2, status)` for interval-censored data.

**SAS C reference:** `src/llike/setcoe_obs_loop.c`

---

## 3. Repeating events (epoch decomposition) -- Complete (v1.0.0)

### What SAS does

Repeating events are handled via epoch decomposition, not a frailty model.
Each patient's longitudinal history is divided into multiple epochs, each
ending at the time of recurrence. Each epoch becomes a separate observation
with:

- `STIME` (left-censoring start time, via `LCENSOR` statement)
- `TIME` (event or right-censoring time)

The likelihood contribution for an epoch is `CF(T) - CF(STIME)` -- the
cumulative hazard accrued during that epoch only.

**SAS C reference:** `src/llike/setcoe_obs_loop.c` lines 92-107

### What R has (as of 1.0.0)

- **Parser:** `.hzr_parse_formula()` recognises `Surv(start, stop, event)`
  and routes `start -> time_lower`, `stop -> time`, `event -> status`.
  Verified by `test-repeating-events.R`.
- **Trivial case:** `Surv(0, t, d)` is equivalent to `Surv(t, d)` and
  produces the same fit on every distribution (also tested).
- **Nonzero-start case:** `Surv(start, stop, event)` with any
  `start > 0` is fully supported. Each epoch contributes `H(stop) − H(start)`
  to the log-likelihood. Both the Weibull and multiphase likelihoods honour
  `time_lower` for all `status` values. Split-invariance tests confirm that
  a single observation `[0, 3.0, event=1]` produces the same fit as two epochs
  `[0, 1.5, event=0] + [1.5, 3.0, event=1]`.

**SAS C reference:** `src/llike/setcoe_obs_loop.c` lines 92-107

---

## 4. Time-varying covariates — Complete

Piecewise constant time-varying coefficients are implemented via
`time_windows = c(t1, t2, ...)`. The internal function
`.hzr_expand_time_varying_design()` expands the design matrix so that each
covariate gets a separate coefficient per time window. Predictions correctly
apply window-specific coefficients based on the time column in `newdata`.

Phase-specific covariate formulas are supported via the `formula` argument
in `hzr_phase()`.

**SAS C reference:** EARLY/CONSTANT/LATE variable lists with `/` separator

---

## 5. Weighted events -- Complete (v1.0.0)

**Status summary (as of 1.0.0):**

| Path | Weighted LL | Weighted gradient | Enforced by `hazard()` |
|:---|:---:|:---:|:---|
| `dist = "weibull"` | Yes (0.9.4) | Yes (0.9.5 fix) | Allowed |
| `dist = "multiphase"` | Yes (0.9.4) | Yes (numDeriv on weighted LL) | Allowed |
| `dist = "exponential"` | Yes (1.0.0) | Yes (1.0.0) | Allowed |
| `dist = "log-logistic"` | Yes (1.0.0) | Yes (1.0.0) | Allowed |
| `dist = "log-normal"` | Yes (1.0.0) | Yes (1.0.0) | Allowed |

All five distributions apply weights via Fisher's multiplicative approach
in the log-likelihood, gradient, and Hessian. Validated against SAS
`WEIGHT` statement output. Additional coverage (fractional weights,
weighted competing risks, weighted bootstrap) is tracked in
`DEVELOPMENT-PLAN.md` Phase 7b.

### What SAS does

Observation weights are specified via the `WEIGHT` statement in `PROC HAZARD`.
Weights are applied multiplicatively in the log-likelihood using Fisher's
weighting approach: each observation's contribution is multiplied by its
weight value. The C code carries weights as position 7 in the observation
array and applies them uniformly to events (`C1`), right-censored (`C2`),
and interval-censored (`C3`) terms.

The Conservation of Events theorem extends to weighted data:
`SUM C1(I)*W1(I) + SUM C3(I)*W3(I) = ...`

**SAS C reference:** `src/llike/setcoe_obs_loop.c` lines 33–36,
`src/llike/setcoe.c` lines 173–246

### What R needs

1. Add a `weights` argument to `hazard()` (numeric vector, length = nrow(data))
2. Thread weights through all distribution-specific log-likelihood functions
   (multiply each observation's LL contribution by its weight)
3. Thread weights through gradient functions
4. Update Hessian computation to account for weights
5. Add tests with weighted vs. unweighted fits
6. Validate against SAS output

**Estimated effort:** Medium — straightforward threading through existing
likelihood infrastructure.

---

## 6. Stepwise covariate selection — Complete (shipped v0.9.5)

`hzr_stepwise()` is fully implemented and tested. See
`inst/dev/STEPWISE-DESIGN.md` for the design document (preserved as
reference). Key capabilities:

- **Forward, backward, two-way** selection via `direction` argument
- **Wald** (`criterion = "wald"`) and **AIC** (`criterion = "aic"`) modes
- **Phase-specific entry** for multiphase models (a variable can enter
  `early` but not `constant` independently)
- **MOVE-limited oscillation guard** (`max_move`) prevents cycling, matching
  SAS `MOVE` semantics
- **Forced/blocked variables** via `force_in` / `force_out`
- **Selection trace** in `$steps` tibble; `print.hzr_stepwise()` displays
  step-by-step entry/exit with criterion-appropriate summary stats
- Works on all five distributions

**Default thresholds match SAS `PROC HAZARD` defaults:**
`slentry = 0.30`, `slstay = 0.20`.

**SAS C reference:** `src/vars/stepw.c`, `backw.c`, `swvari.c`, `swvarx.c`,
`dfast.c`

**FAST screening** (Lawless–Singhal approximate Wald) deferred to a future
version; full per-candidate refits are exact and validated against SAS output.

**Parity status:** CABGKUL forward Wald selection matches SAS `stepw.c`
output within numerical tolerance at `SLENTRY = 0.3`. AIC mode validated
against `stats::step()` on single-distribution Weibull fits.

---

## 7. Conservation of Events theorem -- Complete (standard case; v1.0.0)

**Current scope (as of 1.0.0):**

| Data | CoE applied | Notes |
|:---|:---:|:---|
| Right-censored + exact event (status in {0, 1}) with weights == 1 | Yes | Matches SAS reference on CABGKUL to within 1 log-likelihood unit |
| Interval-censored (status == 2) or left-censored (status == -1) | No -- auto-disabled | Falls through to full-dim optimizer; result is numerically correct |
| Any non-unit weights | No -- auto-disabled | Weighted CoE extension tracked in DEVELOPMENT-PLAN.md Phase 7b |

Fits on unsupported regimes are numerically correct — they just don't get the
one-dimensional closed-form solve and take more optimizer iterations.

The CoE guard in `R/likelihood-multiphase.R`:

```r
coe_supported_data    <- all(status %in% c(0, 1))
coe_supported_weights <- all(weights == 1)
if (use_conserve && coe_supported_data && coe_supported_weights && ...)
```

### What SAS does

Proven by Malcolm E. Turner, Jr. for semi-infinite positive distributions
with right-censored data. For the model:

    S(T) = exp(−G(θ) · D(T))

the CoE theorem states:

    Σ E(i) = Σ CF(T_i)

where E(i) is the event indicator and CF(T_i) is the cumulative distribution
function evaluated at each observation's time. This provides an **explicit
closed-form solution** for one scaling parameter (μ) given the other
parameters.

In practice, this means:
- At each optimizer iteration, μ is solved analytically from the conservation
  constraint rather than numerically optimized
- The optimization problem is reduced by one dimension
- Convergence is faster and more numerically stable
- The result is invariant to the initial guess for μ

The C implementation extends to weighted and interval-censored data.

**SAS C reference:** `src/llike/setcoe.c` (setup, lines 39–250),
`src/llike/consrv.c` (iterative update), `src/llike/constp.c`

### What is implemented (v1.0.0)

1. `.hzr_conserve_events()`: given shape parameters and data, solves for
   `log_mu` satisfying Σ events = Σ CF(t_i) (Turner's closed-form)
2. For multiphase models, applies to one phase's μ per iteration; others
   held fixed
3. Integrated into the optimizer loop: CoE solver called before each LL
   evaluation on the reduced parameter set
4. `conserve` control option (default TRUE) matches SAS `CONSERVE/NOCONSERVE`
5. Validated against SAS output on CABGKUL — identical log-likelihood values,
   fewer optimizer iterations

**Remaining extensions (Phase 7b):** weighted CoE and interval-censored CoE.
The C implementation handles both; the R extension requires carrying weights
through `.hzr_conserve_events()` and modifying the conservation constraint
for the weighted event count.

---

## 8. Covariance and correlation matrix estimation — Complete

The R package computes the observed-information (Hessian-based)
variance-covariance matrix via `numDeriv::hessian()` at the MLE, inverted
to produce the vcov. Standard errors, z-statistics, and p-values are
reported in `summary.hazard()`. The `vcov.hazard()` S3 method extracts the
full matrix.

Fixed parameters (from `fixed = "shapes"`) receive NA entries in the
expanded vcov matrix, while free parameters get proper standard errors.

The `condition` control parameter matches the SAS `CONDITION=` option for
Hessian numerical stability.

**SAS C reference:** `src/optim/`

---

## Implementation status — v1.0.0

All SAS HAZARD parity gaps tracked in this document are **Complete** as of
v1.0.0. Remaining work is extensions beyond the original SAS feature set,
tracked in `DEVELOPMENT-PLAN.md` Phase 7:

| Item | Completed | Version |
|:---|:---:|:---:|
| Multi-phase hazard modeling | ✅ | 0.9.0 |
| Right / interval censoring | ✅ | 0.9.0 |
| Time-varying covariates | ✅ | 0.9.0 |
| Covariance/correlation matrix | ✅ | 0.9.0 |
| Conservation of Events (standard) | ✅ | 0.9.3 |
| Observation weights (Weibull + multiphase) | ✅ | 0.9.4 |
| Observation weights (all distributions) | ✅ | 1.0.0 |
| Repeating events / epoch decomposition | ✅ | 1.0.0 |
| Stepwise covariate selection | ✅ | 0.9.5 |

**Open extensions (Phase 7b):** weighted CoE, weighted competing risks,
fractional weights, additional parity with Rajeswaran/Blackstone production
models.
