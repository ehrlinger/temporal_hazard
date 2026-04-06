# M3 Milestone: Broader Model Support

**Status:** Planning phase. M2 foundation complete; ready to extend distribution family and censoring types.

## Overview

M3 expands the model beyond right-censored Weibull to encompass:
1. Alternative baseline distributions (exponential, log-logistic, log-normal)
2. Interval censoring
3. Time-varying coefficients
4. Comprehensive edge-case parity testing

## Goals

- Support parametric hazard models across multiple distributional families
- Handle diverse censoring mechanisms (right, left, interval)
- Enable time-dependent covariate effects
- Maintain numerical stability and parity with legacy C implementation
- 100% test coverage for new features

## Implementation Strategy

### Phase 1: Distribution Support (Priority 1)

Start with **exponential** (simplest), then **log-logistic** and **log-normal**.

Why this order?
- **Exponential:** No shape parameter; tests baseline architecture without complexity
- **Log-logistic:** Intermediate complexity; commonly used in actuarial work
- **Log-normal:** Widely applied; completes the Trinity of common alternatives

Each distribution requires:
1. Likelihood function `.hzr_logl_<dist>()`
2. Gradient function `.hzr_gradient_<dist>()`
3. Optimizer wrapper `.hzr_optim_<dist>()`
4. `predict()` type support (survival, cumulative_hazard)
5. Unit tests (univariate, multivariate, numerical gradient verification)
6. Golden fixture for parity testing

### Phase 2: Censoring Support (Priority 2)

**Interval censoring** — when event time is known only to lie in [L, U].

Likelihood modification:
- Event observed (δ=1): log h(t|x)
- Right censored (δ=0): log S(t|x)
- Left censored (δ=-1): log[1 - S(t|x)] = log[-log S(t|x)]
- Interval censored (δ=2): log[S(L|x) - S(U|x)]

### Phase 3: Time-Varying Effects (Priority 3)

**Time-varying coefficients** β(t) via stratification.

Approach:
1. Partition follow-up into time windows
2. Fit separate β per window
3. Constrain coefficients across windows (optional shrinkage)

### Phase 4: Parity & Testing (Ongoing)

Expand test harness:
- Boundary conditions (t → 0, t → ∞)
- Convergence failures and recoveries
- High-dimensional edge cases (p > n)
- Censoring mechanism robustness

---

## Current Implementation Status

### Completed (from M2)

- ✅ Weibull likelihood and gradient
- ✅ L-BFGS-B optimizer with bounds
- ✅ survival and cumulative_hazard predictions
- ✅ Golden fixtures and parity tests (Weibull only)
- ✅ 83 passing test suite

### In Progress

(None yet; awaiting kickoff)

### Pending

- [ ] Exponential distribution core
- [ ] Log-logistic distribution core
- [ ] Log-normal distribution core
- [ ] Interval censoring
- [ ] Time-varying coefficients
- [ ] Extended parity suite

---

## Distribution Mathematical Summary

### Exponential: h(t) = λ

```
h(t) = λ
S(t) = exp(-λt)
H(t) = λt
ℓ = sum(δ_i log λ) - λ sum(t_i)
```

**Parameters:** λ > 0

### Log-Logistic: h(t) = (α·β·t^(β-1)) / (1 + α·t^β)

```
S(t) = 1 / (1 + α·t^β)
H(t) = log(1 + α·t^β)
```

**Parameters:** α > 0 (scale), β > 0 (shape)

### Log-Normal: h(t) = φ(log(t) | μ, σ) / [σ · t · (1 - Φ(log(t) | μ, σ))]

```
S(t) = 1 - Φ((log(t) - μ) / σ)
```

**Parameters:** μ ∈ ℝ (location), σ > 0 (scale)

---

## Next Action

Start with **Task 1: Create distribution module structure and exponential core**.

Expected output:
- `R/likelihood-exponential.R` — likelihood, gradient, optimizer
- `tests/testthat/test-exponential-dist.R` — 10+ unit tests
- Update `hazard()` to dispatch on dist = "exponential"
- Golden fixture: `inst/fixtures/hz_exponential.rds`
