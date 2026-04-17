# SAS HAZARD Parity Gap Analysis

**Date:** 2026-04-16
**Package version:** 0.9.4
**Branch:** feature/weights-repeating-events

This document tracks feature parity between the original C/SAS HAZARD system
(Blackstone, Naftel, Turner 1986; maintained by Rajeswaran & Ehrlinger) and
the TemporalHazard R package.

## Summary

| Capability | R Status | SAS Reference |
|:---|:---:|:---|
| Multi-phase hazard modeling | Complete | `src/llike/`, `src/decomp/` |
| Right and interval censoring | Complete | `src/llike/setcoe_obs_loop.c` |
| Repeating events (epoch decomposition) | Planned | `src/llike/setcoe_obs_loop.c` (STIME) -- parser accepts `Surv(start, stop, event)` but LL ignores `start` for counting-process rows (see sec. 3) |
| Time-varying covariates | Complete | EARLY/CONSTANT/LATE variable lists |
| Weighted events | Partial | `src/llike/setcoe_obs_loop.c` -- Weibull + multiphase only; exp / log-logistic / log-normal pending (see sec. 5) |
| Stepwise covariate selection | Planned | `src/vars/stepw.c`, `backw.c`, `swvari.c` |
| Conservation of Events theorem | Complete | `src/llike/setcoe.c`, `consrv.c` |
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

## 3. Repeating events (epoch decomposition) -- Planned (narrowed in v0.9.5)

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

### What R has (as of 0.9.5)

- **Parser:** `.hzr_parse_formula()` recognises `Surv(start, stop, event)`
  and routes `start -> time_lower`, `stop -> time`, `event -> status`.
  Verified by `test-repeating-events.R`.
- **Trivial case:** `Surv(0, t, d)` is equivalent to `Surv(t, d)` and
  produces the same fit on every distribution (also tested).
- **Nonzero-start case:** `Surv(start, stop, event)` with any
  `start > 0` is **rejected at `hazard()`** in v0.9.5. The 0.9.4 NEWS
  claimed epochs contributed `H(stop) - H(start)` to the likelihood,
  but the downstream likelihoods (Weibull single-dist and multiphase)
  only honour `time_lower` for interval-censored rows (`status == 2`).
  Counting-process rows (`status` in {0, 1}) were silently scored as
  `H(stop)` alone -- an observation at `[1.5, 3.0, event=1]` got the
  same LL contribution as `[0, 3.0, event=1]`. The sub-case is
  blocked with an explicit error until the wire-up is finished.

### What is still needed (Phase 4f)

1. In the Weibull LL (`.hzr_logl_weibull`): swap
   `cumhaz_event[idx_event] -> cumhaz_event[idx_event] - cumhaz_lower[idx_event]`
   in the event term and
   `cumhaz_event[idx_right] -> cumhaz_event[idx_right] - cumhaz_lower[idx_right]`
   in the right-censored term, honouring `time_lower` for all
   `status` values, not just `status == 2`.
2. Same swap in `.hzr_logl_multiphase()`.
3. Analytic gradient updates (straightforward: all terms pick up an
   additional `-cumhaz_lower` contribution).
4. Remove the `hazard()` guard on counting-process data and re-enable
   the split-invariance tests in `test-repeating-events.R`.
5. Add a vignette section on epoch-decomposed longitudinal data and a
   parity test against a SAS reference fit.

**Estimated effort:** Medium. The likelihood math is already in place
for the interval-censored path; this is a small surgical extension plus
analytic gradient bookkeeping plus test + doc work.

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

## 5. Weighted events -- Partial

**Status summary (as of 0.9.5):**

| Path | Weighted LL | Weighted gradient | Enforced by `hazard()` |
|:---|:---:|:---:|:---|
| `dist = "weibull"` | Yes (0.9.4) | Yes (0.9.5 fix) | Allowed |
| `dist = "multiphase"` | Yes (0.9.4) | Yes (numDeriv on weighted LL) | Allowed |
| `dist = "exponential"` | **No** | **No** | `hazard()` errors |
| `dist = "log-logistic"` | **No** | **No** | `hazard()` errors |
| `dist = "log-normal"` | **No** | **No** | `hazard()` errors |

The three unfinished single-distribution paths accept `weights` as a
formal argument for signature consistency but never apply it inside the
likelihood. Rather than return a silently-unweighted MLE, `hazard()`
raises an error when `weights` is supplied with one of those
distributions. Full wire-up is tracked in `DEVELOPMENT-PLAN.md` Phase 8.

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

## 6. Stepwise covariate selection — Planned

### What SAS does

The C/SAS HAZARD system has a full stepwise variable selection engine:

- **Forward selection:** Add variables sequentially; remove those failing
  SLSTAY criterion
- **Backward elimination:** Remove variables; re-examine for additions via
  SLENTRY
- **Two-way stepwise:** Alternate forward and backward phases
- **FAST screening:** Approximate Wald test updates using the Lawless &
  Singhal method for efficient initial screening

Key parameters:
- `SLENTRY` (default 0.3): p-value threshold for entry
- `SLSTAY` (default 0.2 forward, 0.05 backward): p-value threshold for
  retention
- `MOVE`: Maximum times a variable can oscillate in/out (prevents cycling)
- `ORDER`: Forced entry sequence
- `MAXSTEPS`, `MAXVARS`: Iteration and model size limits

Selection is performed **within** the hazard model estimation, not as a
preprocessing step. This is important because it accounts for the multiphase
structure when evaluating covariate significance.

**SAS C reference:** `src/vars/stepw.c`, `backw.c`, `swvari.c`, `swvarx.c`,
`dfast.c`

### What R needs

1. Wald test infrastructure for individual covariate effects within the
   multiphase framework (z-statistics from the vcov matrix)
2. Forward step: fit model with each candidate variable added, select the
   one with the smallest p-value below SLENTRY
3. Backward step: test each included variable, remove the one with the
   largest p-value above SLSTAY
4. Stepping loop with MOVE limits and convergence detection
5. FAST screening mode using approximate Wald updates (optional optimization)
6. Phase-specific variable selection (variables can enter different phases
   independently)
7. User-facing API: `hzr_stepwise()` or `step()` method for hazard objects
8. Printing/reporting of selection trace

**Estimated effort:** Large — this is the most complex missing feature.
The basic forward/backward logic is straightforward, but the multiphase
phase-specific selection and FAST screening are substantial.

---

## 7. Conservation of Events theorem — Complete (v0.9.3)

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

### What R needs

1. Implement the CoE constraint as a function: given shape parameters and
   data, solve for the log_mu that satisfies Σ events = Σ CF(t_i)
2. For multiphase models, this applies to one phase's μ at a time (the
   others are held fixed or solved iteratively)
3. Integrate into the optimizer loop: at each iteration, call the CoE solver
   before evaluating the likelihood on the reduced parameter set
4. Handle the extension to interval-censored and weighted data
5. Add a `conserve` control option (default TRUE) matching SAS
   `CONSERVE/NOCONSERVE`
6. Validate against SAS output — CoE should yield identical log-likelihood
   values with fewer iterations

**Estimated effort:** Medium-large — the math is well-defined (Turner's
proof), but integrating it into the multiphase optimizer loop requires
careful handling of the reduced-dimensional optimization and the interaction
with fixed shape parameters.

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

## Implementation priority

1. **Stepwise selection** — critical for exploratory clinical research
   workflows; large effort; high user-facing impact

Completed:
- Conservation of Events (v0.9.3) — integrated into multiphase optimizer
- Observation weights (v0.9.4) — Fisher weighting threaded through likelihood,
  gradient, and Hessian
- Repeating events (v0.9.4) — `Surv(start, stop, event)` epoch decomposition
