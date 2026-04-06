test_that("Golden fixture: hz.death.AVC univariable Weibull", {
  skip("Parity harness under development (M2)")
  
  # This test will:
  # 1. Load AVC death data from examples/data/avc
  # 2. Run C hazard binary with same parameters as SAS analysis
  # 3. Fit R hazard() model with same inputs
  # 4. Compare key outputs (log-likelihood, parameter estimates, SE, etc.)
  # 5. Assert parity within tolerance
  
  # Data loading stub
  # avc <- read.table(...)
  
  # Golden fixture generation
  # golden <- .hzr_generate_golden_fixture(...)
  
  # R model fitting
  # fit_r <- hazard(...)
  
  # Parity assertions
  # expect_equal(fit_r$fit$objective, golden$parsed$log_lik, tolerance = 1e-4)
})

test_that("Golden fixture: hm.death.AVC multivariable Weibull", {
  skip("Parity harness under development (M2)")
  
  # Multivariable model with covariates
})

test_that("Parity tolerance for numerical differences", {
  skip("Parity harness under development (M2)")
  
  # Define acceptable tolerance ranges for:
  # - Log-likelihood (absolute)
  # - Parameter estimates (relative %)
  # - Standard errors (relative %)
  # - Gradient norm at convergence
})
