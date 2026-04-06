# M2 Milestone: Model Core & Likelihood Evaluation

**Status:** Early implementation; framework in place; ready for gradient porting and golden fixture generation.

## What's Complete

### 1. Weibull Likelihood Foundation
- **File:** `R/likelihood-weibull.R`
- **Core Functions:**
  - `.hzr_logl_weibull()` — Log-likelihood evaluation with numerical stability safeguards
  - `.hzr_gradient_weibull()` — Placeholder for score vector (ready for C port)
  - `.hzr_optim_weibull()` — L-BFGS-B optimizer with box constraints (mu > 0, nu > 0)

### 2. API Extensions
- **`hazard()` now accepts `fit = TRUE/FALSE`** (default FALSE for M1 compatibility)
  - When `fit=TRUE` and `theta` provided: triggers ML estimation via optimizer
  - Returns fitted `$fit$theta`, `$fit$objective` (log-lik), `$fit$converged` flags
  
- **New S3 Methods:**
  - `coef.hazard()` — Extract estimated coefficients
  - `vcov.hazard()` — Extract variance-covariance matrix (for SE computation)

### 3. C Binary Integration
- **Compiled Reference:** `inst/bin/hazard` (937 KB, freshly built from source)
- **Purpose:** Generate golden fixtures for parity validation
- **Helpers (`R/parity-helpers.R`):**
  - `.hzr_get_hazard_binary()` — Locate and validate binary
  - `.hzr_run_hazard_binary()` — Execute with data; parse output (scaffold)
  - `.hzr_generate_golden_fixture()` — Save reference outputs as test fixtures
  - `.hzr_parse_hazard_output()` — Placeholder for output parsing

### 4. Parity Test Harness
- **File:** `tests/testthat/test-parity-core.R`
- **Status:** Skipped stubs pending fixture implementation
- **Plan:** 
  1. Load AVC and other real datasets
  2. Run C binary on same inputs as R fit
  3. Assert parameter estimates, log-lik, etc. within tolerance
  4. Define tolerance table (% diffs for params, absolute for log-lik)

## What's In Progress

### Gradient Calculation
The `.hzr_gradient_weibull()` function is currently a stub. To complete:
1. Compute partial derivatives of likelihood w.r.t. theta (shape and covariate parameters)
2. Port math from C implementation (`src/llike/hzd_calc_parm_drv.c`, etc.)
3. Return vector matching theta parametrization
4. Optional: Implement numerical gradient as fallback

### Optimizer Refinement
- Current: L-BFGS-B with basic penalty for infeasible params
- Next: Enhanced convergence diagnostics (gradient norm, relative param change)
- Future: Consider trust-region methods for better numerical stability

## What's Pending (M2 Phase 2)

### Golden Fixture Generation
1. **Univariable Models** (shape estimation only)
   - `hz.death.AVC` — Weibull fit to AVC death data
   - `hz.deadp.KUL` — Weibull fit to KUL dataset
   
2. **Multivariable Models** (shape + covariates)
   - `hm.death.AVC` — Full multivariable with EARLY/CONSTANT phases
   - `hm.deadp.VALVES` — Valve surgery dataset

### predict() Type Expansion
Current: `"linear_predictor"`, `"hazard"`
Needed for M2:
- `"survival"` — S(t | x) = exp(-cumulative hazard)
- `"cumulative_hazard"` — Cumulative hazard at time t

### Parity Validation
Define tolerance table:
```
Parameter Type       | Tolerance
---------------------|----------
Log-likelihood       | ±1e-4 (absolute)
Shape parameters     | ±0.1% (relative)
Covariate coefs      | ±0.2% (relative)
Standard errors      | ±0.5% (relative)
Gradient norm at sol | < 1e-5
```

## Quick Integration Example

```r
library(TemporalHazard)

# Load data
avc <- read.table("examples/data/avc", header = FALSE)

# Fit univariable Weibull (no covariates)
fit_uni <- hazard(
  time = avc$INT_DEAD,
  status = avc$DEAD,
  x = NULL,
  theta = c(mu = 0.3, nu = 1.4),  # Starting values
  dist = "weibull",
  fit = TRUE,                      # Enable optimization
  control = list(maxit = 1000, reltol = 1e-5)
)

# Extract fitted values
coef(fit_uni)                      # Estimated [mu, nu]
vcov(fit_uni)                      # Variance-covariance matrix

# Predict on original data
eta <- predict(fit_uni, type = "linear_predictor")  # [0] (no covariates)
haz <- predict(fit_uni, type = "hazard")            # [1] (baseline)
```

## Files Updated

- `R/likelihood-weibull.R` — NEW: Core likelihood and optimizer
- `R/hazard_api.R` — Updated: Added `fit` parameter, coef/vcov methods
- `R/parity-helpers.R` — NEW: Golden fixture infrastructure
- `tests/testthat/test-parity-core.R` — NEW: Parity validation stubs
- `inst/bin/hazard` — NEW: Compiled C reference binary
- `NAMESPACE` — Updated: Exports `coef.hazard`, `vcov.hazard`

## Known Issues & To-Do

- [ ] Gradient calculation needs full port from C code
- [ ] Numerical Hessian computation could be faster (consider analytical form)
- [ ] Golden fixture XPORT format parsing not yet implemented
- [ ] Multi-phase model (EARLY/CONSTANT/LATE) structure not yet formalized
  
## Testing Status

```
Status: 0 ERRORS, 4 WARNINGS, 2 NOTES
Tests: 27 passed, 3 skipped (M2 parity harness)
Vignettes: 9 converted to Quarto (no executable code yet — awaiting M2 examples)
```

## Next Immediate Steps

1. **Implement gradient** (1-2 hours)
   - Symbolic or numerical differentiation of likelihood
   - Validate against C code gradients
   
2. **Test fitting on small synthetic data** (30 min)
   - Confirm convergence behavior
   - Tune L-BFGS-B tolerances
   
3. **Generate first golden fixture** (1 hour)
   - Run C binary on AVC univariable
   - Capture output; store as `.rds`
   - Implement `.hzr_parse_hazard_output()` for real output format
   
4. **Implement parity tests** (2 hours)
   - Load fixture and fitted R model
   - Assert parity within tolerance
   - Document any discrepancies

---

**Est. Full M2 Completion:** 4-6 additional hours of focused porting and fixture generation.
