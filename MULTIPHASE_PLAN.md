# Multiphase Hazard Model — Implementation Plan

Status: DRAFT — ready for review before implementation begins.

## Design Principles

1. **Generalized N-phase** — every phase uses the same `decompos()` engine.
   No hard-coded 3-phase limit.
2. **Prerelease** — the existing single-distribution API can change freely.
   No backward-compatibility shims.
3. **Modern R** — S3 classes, roxygen2, formula interface, clean NAMESPACE.
4. **SAS/C bridge** — parameter names recognizable to HAZARD veterans.
   `hzr_argument_mapping()` extended to cover multiphase parameters.

## Mathematical Foundation

### Generalized Temporal Decomposition (from mixhazard)

Every phase uses the parametric family `decompos(t; t_half, nu, m)`, which
returns three quantities:

- **G(t)** — cumulative distribution (CDF)
- **g(t)** — density: the "early" temporal pattern
- **h(t) = g(t) / (1 − G(t))** — hazard: the "late" temporal pattern

The three parameters control the shape:

- **t_half** — half-life (time at which G = 0.5); must be > 0
- **nu** (η in mixhazard) — exponent of time; controls rate dynamics
- **m** (γ in mixhazard) — shape exponent; controls distributional form

Six valid sign combinations span the shape space:

| Case | m   | nu  | Behavior              |
|------|-----|-----|-----------------------|
| 1    | > 0 | > 0 | Standard sigmoidal    |
| 1L   | = 0 | > 0 | Exponential-like      |
| 2    | < 0 | > 0 | Heavy-tailed          |
| 2L   | < 0 | = 0 | Exponential decay     |
| 3    | > 0 | < 0 | Bounded cumulative    |
| 3L   | = 0 | < 0 | Bounded exponential   |

Invalid: both m < 0 AND nu < 0.

### SAS/C Parameter Name Mapping

The original C code uses different parameter names per phase. The
generalized decomposition unifies them:

| C/SAS Phase | C Parameters            | decompos Equivalent      |
|-------------|-------------------------|--------------------------|
| Early (G1)  | DELTA, RHO/THALF, NU, M | t_half, nu, m            |
| Constant    | (none — G₂(t) = t)     | type = "constant"        |
| Late (G3)   | TAU, GAMMA, ALPHA, ETA  | t_half, nu, m            |

The 4-parameter C/SAS parameterizations (DELTA/RHO/NU/M for early,
TAU/GAMMA/ALPHA/ETA for late) collapse onto the 3-parameter decompos
family. The C DELTA parameter controlled a time transformation
`B(t) = (exp(δt) − 1)/δ` that decompos absorbs into the shape.

### Additive Multiphase Model

Total cumulative hazard decomposes additively across J phases:

```
H(t | x) = Σ_{j=1}^{J}  μ_j(x) · Φ_j(t; t_half_j, nu_j, m_j)
```

Where:

- `μ_j(x) = exp(α_j + x · β_j)` — phase-specific log-linear scaling
  with covariates (SAS: `MU_j` and `B_j` coefficients)
- `Φ_j(t)` — temporal shape for phase j:
  - `"cdf"` → Φ = G(t). Bounded [0, 1]. Early risk that resolves.
  - `"hazard"` → Φ = ∫₀ᵗ h(s) ds. Monotone increasing. Late/aging risk.
  - `"constant"` → Φ = t. Flat background rate. (SAS: G2 phase.)

Total hazard:

```
h(t | x) = Σ_{j=1}^{J}  μ_j(x) · φ_j(t)
```

where φ_j = dΦ_j/dt.

Survival:

```
S(t | x) = exp(−H(t | x))
```

## API Design

### Phase Specification

```r
#' Specify a hazard phase
#'
#' @param type Phase type: "cdf" (early), "hazard" (late), or "constant".
#' @param t_half Initial half-life (> 0). Ignored for "constant".
#' @param nu Initial time exponent. Ignored for "constant".
#'   SAS early: NU; SAS late: maps from GAMMA/ALPHA/ETA.
#' @param m Initial shape. Ignored for "constant".
#'   SAS early: M; SAS late: maps from GAMMA/ALPHA/ETA.
#' @param formula One-sided formula for phase-specific covariates,
#'   e.g. ~ age + nyha. NULL = same covariates as global formula.
#' @return An S3 "hzr_phase" object.
#' @export
hzr_phase <- function(type = c("cdf", "hazard", "constant"),
                      t_half = 1, nu = 1, m = 0,
                      formula = NULL)
```

### hazard() Interface

```r
# ── Single-phase (standard parametric) ──────────────────────────────
fit <- hazard(
  Surv(time, status) ~ age + nyha,
  data  = dat,
  dist  = "weibull",
  theta = c(mu = 0.25, nu = 1.1, beta1 = 0, beta2 = 0),
  fit   = TRUE
)

# ── Two-phase decomposition ─────────────────────────────────────────
fit2 <- hazard(
  Surv(time, status) ~ age + nyha,
  data   = dat,
  dist   = "multiphase",
  phases = list(
    early = hzr_phase("cdf",    t_half = 0.5, nu = 2, m = 0),
    late  = hzr_phase("hazard", t_half = 5,   nu = 1, m = 0)
  ),
  fit = TRUE
)

# ── Three-phase (classic HAZARD pattern) ────────────────────────────
fit3 <- hazard(
  Surv(time, status) ~ age,
  data   = dat,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf",      t_half = 0.5, nu = 2, m = 0),
    constant = hzr_phase("constant"),
    late     = hzr_phase("hazard",   t_half = 5,   nu = 1, m = 0)
  ),
  fit = TRUE
)

# ── Phase-specific covariates ───────────────────────────────────────
fit4 <- hazard(
  Surv(time, status) ~ age + nyha + shock,
  data   = dat,
  dist   = "multiphase",
  phases = list(
    early = hzr_phase("cdf",    t_half = 0.5, nu = 2, m = 0,
                      formula = ~ age + shock),
    late  = hzr_phase("hazard", t_half = 5,   nu = 1, m = 0,
                      formula = ~ age + nyha)
  ),
  fit = TRUE
)
```

Named phases (`early = ...`, `late = ...`) are used in output labels.
Unnamed phases get auto-labeled `phase_1`, `phase_2`, etc.

### predict() with Decomposition

```r
# Total survival / cumulative hazard (same API as single-phase)
predict(fit2, newdata = grid, type = "survival")
predict(fit2, newdata = grid, type = "cumulative_hazard")

# Per-phase decomposition — returns data frame with phase columns
predict(fit2, newdata = grid, type = "cumulative_hazard",
        decompose = TRUE)
# Returns: time, total, early, late  (one column per named phase)
```

### summary() Output

```
Multiphase hazard model (2 phases)

Phase 1: early (cdf)
  Shape:  t_half = 0.482 (0.031), nu = 2.14 (0.18), m = 0.12 (0.05)
  Scale:  alpha = -1.23 (0.41)
  Covariates:
    age   0.032 (0.008)  ***
    shock 0.891 (0.214)  ***

Phase 2: late (hazard)
  Shape:  t_half = 5.31 (0.72), nu = 0.98 (0.11), m = 0.01 (0.03)
  Scale:  alpha = -3.45 (0.33)
  Covariates:
    age   0.041 (0.006)  ***
    nyha  0.156 (0.072)  *

Log-likelihood: -287.3 on 12 parameters
Convergence: 0 (success)
```

## Parameter Vector Layout

### Internal (estimation scale)

For each phase j, the parameter block is:

```
[log_mu_j, log_t_half_j, nu_j, m_j, beta_j_1, ..., beta_j_pj]
```

- `log_mu_j = α_j` — log-scale intercept (unconstrained)
- `log_t_half_j` — log-transformed half-life (ensures > 0)
- `nu_j`, `m_j` — unconstrained shape parameters
- `beta_j_k` — covariate coefficients (unconstrained)
- "constant" phases omit `log_t_half`, `nu`, `m`

### User-facing (reporting scale)

Back-transformed for interpretability:

- `mu_j = exp(α_j)` — baseline cumulative hazard scaling
- `t_half_j = exp(log_t_half_j)` — half-life in original time units
- `nu_j`, `m_j` — reported as-is
- `beta_j_k` — reported as-is (log-hazard scale)

## Implementation Steps

### Step 1: Core Decomposition Engine

**New file: `R/decomposition.R`**

Port `decompos()` from mixhazard with these changes:

- Rename to `hzr_decompos()` with `@export`
- Use parameter names `t_half`, `nu`, `m` (matching SAS early-phase naming)
- Return named list: `G` (CDF), `g` (density), `h` (hazard)
- Add vectorized `hzr_phase_cumhaz(time, t_half, nu, m, type)` — returns
  Φ(t) for the given phase type
- Add vectorized `hzr_phase_hazard(time, t_half, nu, m, type)` — returns
  φ(t) = dΦ/dt
- Numerical guards: clamp time > .Machine$double.xmin, handle G(t) → 1

**New file: `tests/testthat/test-decomposition.R`**

- Validate all 6 cases against mixhazard::decompos() output
- Test boundary cases (nu = 0, m = 0)
- Test invalid case (m < 0 && nu < 0) throws error
- Test monotonicity of G(t)
- Test h(t) = g(t) / (1 − G(t)) identity

### Step 2: Phase Specification

**New file: `R/phase-spec.R`**

- `hzr_phase()` constructor — validates type, stores initial values
- `print.hzr_phase()` — compact display
- `is_hzr_phase()` — type check
- `.hzr_phase_n_shape()` — returns 3 for cdf/hazard, 0 for constant
- `.hzr_phase_theta_names()` — generates named vector labels

### Step 3: Multiphase Likelihood

**New file: `R/likelihood-multiphase.R`**

Following the existing likelihood file pattern:

- `.hzr_logl_multiphase(theta, time, status, ..., phases)` — negative
  log-likelihood using the additive model
- `.hzr_gradient_multiphase(theta, time, status, ..., phases)` — analytical
  score vector (sum of per-phase contributions)
- `.hzr_optim_multiphase(time, status, ..., x, theta_start, control, phases)`
  — multi-start L-BFGS-B optimizer

Optimization details:
- Internal parameterization: `(log_mu, log_t_half, nu, m, beta)`
- L-BFGS-B bounds: log_t_half > log(1e-4), others unconstrained
- Feasibility guard: reject theta where both m < 0 and nu < 0
- Multi-start: `control$n_starts` (default 5) random perturbations around
  user-supplied starting values
- Delta method for vcov back-transformation

### Step 4: Wire into hazard() API

**Modify: `R/hazard_api.R`**

- Add `phases` parameter to `hazard()`
- When `dist = "multiphase"`: validate phases, build theta_start from phase
  specs, dispatch to `.hzr_optim_multiphase()`
- Store `phases` in `object$spec$phases`
- Update `print.hazard()`, `summary.hazard()` for multiphase display

### Step 5: Extend predict() and summary()

**Modify: `R/hazard_api.R`**

- `predict.hazard()`: add `decompose = FALSE` parameter. When TRUE and
  dist = "multiphase", return a data frame with per-phase cumulative hazard
  columns alongside the total.
- `summary.hazard()`: group parameters by phase with interpretable labels.

### Step 6: Update argument_mapping.R

**Modify: `R/argument_mapping.R`**

Add rows for multiphase parameters:

| SAS/C Input | R Parameter | Notes |
|-------------|-------------|-------|
| MU_1, MU_2, MU_3 | Phase-specific mu (via hzr_phase) | exp(alpha_j) |
| THALF | t_half | Half-life in hzr_phase() |
| NU (early) | nu | Time exponent in hzr_phase() |
| M (early) | m | Shape in hzr_phase() |
| TAU (late) | t_half | Maps to half-life concept |
| GAMMA (late) | nu, m | Collapsed into decompos params |
| ALPHA (late) | m | Part of shape mapping |
| ETA (late) | nu | Part of exponent mapping |
| DELTA (early) | (absorbed) | Time transform absorbed by decompos |

### Step 7: Phase Decomposition Plot

Add to hazard() Rd example and create a new vignette section showing the
multiphase decomposition — total curve with per-phase components, using
hvtiPlotR with color-coded phases and a legend.

### Step 8: Tests

| Test File | Coverage |
|-----------|----------|
| test-decomposition.R | hzr_decompos() all 6 cases, parity with mixhazard |
| test-phase-spec.R | hzr_phase() construction, validation, theta names |
| test-multiphase-likelihood.R | Known-parameter evaluation, gradient vs numerical |
| test-multiphase-api.R | End-to-end: 1, 2, 3 phases; predict; summary |
| test-multiphase-parity.R | 3-phase fit matches C code output on test data |

## File Summary

| File | Action | Description |
|------|--------|-------------|
| R/decomposition.R | NEW | hzr_decompos(), phase cumhaz/hazard |
| R/phase-spec.R | NEW | hzr_phase() constructor and helpers |
| R/likelihood-multiphase.R | NEW | Likelihood, gradient, optimizer |
| R/hazard_api.R | MODIFY | phases param, multiphase dispatch, predict decompose |
| R/argument_mapping.R | MODIFY | Add multiphase SAS/C mappings |
| tests/testthat/test-decomposition.R | NEW | Decomposition unit tests |
| tests/testthat/test-phase-spec.R | NEW | Phase spec tests |
| tests/testthat/test-multiphase-likelihood.R | NEW | Likelihood tests |
| tests/testthat/test-multiphase-api.R | NEW | API integration tests |
| tests/testthat/test-multiphase-parity.R | NEW | C code parity tests |
| vignettes/multiphase-example.qmd | NEW | Multiphase decomposition vignette |

## Build Order

Each step is independently testable:

1. `decomposition.R` + test-decomposition.R
2. `phase-spec.R` + test-phase-spec.R
3. `likelihood-multiphase.R` + test-multiphase-likelihood.R
4. `hazard_api.R` changes + test-multiphase-api.R
5. predict/summary extensions
6. argument_mapping.R updates
7. Vignette and Rd example plots
8. test-multiphase-parity.R (requires C binary or reference output)
