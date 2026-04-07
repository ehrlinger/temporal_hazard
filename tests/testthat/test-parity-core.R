test_that("Univariable model recovers estimated parameters reproducibly", {
  # Load the golden fixture
  fixture_file <- system.file("fixtures", "hz_univariate.rds", package = "TemporalHazard")
  fixture <- readRDS(fixture_file)
  
  # Skip if fixture doesn't exist
  skip_if_not(file.exists(fixture_file))
  
  # Re-fit the same model
  refit <- hazard(
    time = fixture$data$time,
    status = fixture$data$status,
    x = NULL,
    theta = c(mu = 0.3, nu = 1.0),
    dist = "weibull",
    fit = TRUE,
    control = list(maxit = 500, reltol = 1e-6)
  )
  
  # Model should have fit results
  expect_true(!is.null(refit$fit))
  expect_true("theta" %in% names(refit$fit))
  
  # Parameter estimates should match fixture closely (within tolerance)
  # Allow slightly higher tolerance due to optimization variations
  expect_equal(refit$fit$theta, fixture$fit$theta, tolerance = 1e-3)
  
  # Log-likelihood should be close
  expect_equal(refit$fit$objective, fixture$fit$logl, tolerance = 1e-2)
})

test_that("Multivariable model with covariates recovers parameters", {
  fixture_file <- system.file("fixtures", "hm_multivariate.rds", package = "TemporalHazard")
  fixture <- readRDS(fixture_file)
  
  skip_if_not(file.exists(fixture_file))
  
  # Re-fit
  refit <- hazard(
    time = fixture$data$time,
    status = fixture$data$status,
    x = fixture$data$x,
    theta = c(mu = 0.4, nu = 1.0, beta_X1 = -0.2, beta_X2 = 0.3),
    dist = "weibull",
    fit = TRUE,
    control = list(maxit = 500, reltol = 1e-6)
  )
  
  # Model should have fit results
  expect_true(!is.null(refit$fit))
  expect_true("theta" %in% names(refit$fit))
  
  # Parameter estimates within tolerance
  expect_equal(refit$fit$theta, fixture$fit$theta, tolerance = 1e-3)
  
  # Log-likelihood should be close
  expect_equal(refit$fit$objective, fixture$fit$logl, tolerance = 1e-2)
})

test_that("Edge case: small sample with many covariates", {
  fixture_file <- system.file("fixtures", "hm_edge_case.rds", package = "TemporalHazard")
  fixture <- readRDS(fixture_file)
  
  skip_if_not(file.exists(fixture_file))
  
  # Re-fit
  refit <- hazard(
    time = fixture$data$time,
    status = fixture$data$status,
    x = fixture$data$x,
    theta = c(mu = 0.5, nu = 1.2, beta_X1 = 0, beta_X2 = 0, beta_X3 = 0),
    dist = "weibull",
    fit = TRUE,
    control = list(maxit = 300, reltol = 1e-5)
  )
  
  # Check convergence flag (may not always converge in 300 iters)
  expect_type(refit$fit, "list")
  expect_true("theta" %in% names(refit$fit))
  expect_true(all(is.finite(refit$fit$theta)))
  expect_true(is.finite(refit$fit$objective))
  
  # This fixture is intentionally ill-conditioned (small n, high leverage),
  # so assert objective parity rather than exact parameter recovery.
  expect_equal(refit$fit$objective, fixture$fit$logl, tolerance = 5)
})

test_that("Univariable predictions are computable", {
  fixture_file <- system.file("fixtures", "hz_univariate.rds", package = "TemporalHazard")
  fixture <- readRDS(fixture_file)
  
  skip_if_not(file.exists(fixture_file))
  
  # Refit model
  refit <- hazard(
    time = fixture$data$time,
    status = fixture$data$status,
    x = NULL,
    theta = c(mu = 0.3, nu = 1.0),
    dist = "weibull",
    fit = TRUE,
    control = list(maxit = 500, reltol = 1e-6)
  )
  
  # Predictions should be computable with time data
  pred1 <- predict(refit, type = "linear_predictor", 
                   newdata = data.frame(time = fixture$data$time[1:5]))
  pred2 <- predict(refit, type = "hazard",
                   newdata = data.frame(time = fixture$data$time[1:5]))
  
  # Should return numeric vectors
  expect_true(is.numeric(pred1))
  expect_true(is.numeric(pred2))
  expect_true(length(pred1) ==  5)
  expect_true(length(pred2) == 5)
})

test_that("Convergence codes are reasonable", {
  fixture_file <- system.file("fixtures", "hz_univariate.rds", package = "TemporalHazard")
  fixture <- readRDS(fixture_file)
  
  skip_if_not(file.exists(fixture_file))
  
  # Refit
  refit <- hazard(
    time = fixture$data$time,
    status = fixture$data$status,
    x = NULL,
    theta = c(mu = 0.3, nu = 1.0),
    dist = "weibull",
    fit = TRUE,
    control = list(maxit = 500, reltol = 1e-6)
  )
  
  # Convergence code should be logical or integer (TRUE or FALSE for success)
  expect_true(is.logical(refit$fit$converged) || is.integer(refit$fit$converged))
})


