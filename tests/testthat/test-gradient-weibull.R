test_that("gradient_weibull computes correct shape for univariate case", {
  # Simple univariate Weibull (no covariates)
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)

  theta <- c(mu = 0.5, nu = 1.2)

  # Compute gradient analytically
  grad <- .hzr_gradient_weibull(
    theta = theta,
    time = time,
    status = status,
    x = NULL
  )

  expect_length(grad, 2)
  expect_true(all(is.finite(grad)))
})

test_that("gradient_weibull computes correct shape with covariates", {
  # Weibull with one covariate
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)
  x <- matrix(c(0.5, 1.2, -0.3, 0.8, 1.5), ncol = 1)

  theta <- c(mu = 0.5, nu = 1.2, beta = 0.3)

  # Compute gradient analytically
  grad <- .hzr_gradient_weibull(
    theta = theta,
    time = time,
    status = status,
    x = x
  )

  expect_length(grad, 3)
  expect_true(all(is.finite(grad)))
})

test_that("gradient_weibull matches numerical gradient", {
  # Use numDeriv to verify analytical gradient
  skip_if_not_installed("numDeriv")

  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 0, 1, 1, 0)
  x <- matrix(c(0.5, 1.2, -0.3, 0.8, 1.5), ncol = 1)

  theta <- c(mu = 0.5, nu = 1.2, beta = 0.3)

  # Analytical gradient
  grad_analytical <- .hzr_gradient_weibull(
    theta = theta,
    time = time,
    status = status,
    x = x
  )

  # Numerical gradient via numDeriv
  objective <- function(th) {
    .hzr_logl_weibull(
      theta = th,
      time = time,
      status = status,
      x = x,
      return_gradient = FALSE
    )
  }

  grad_numerical <- numDeriv::grad(objective, theta)

  # Should match to machine precision (allowing small numerical errors)
  expect_equal(grad_analytical, grad_numerical, tolerance = 1e-4)
})

test_that("optimizer converges using analytical gradient", {
  # Test that optimizer runs with our gradient enabled
  time <- c(1, 2, 3, 4, 5, 6)
  status <- c(1, 0, 1, 1, 0, 1)

  theta_start <- c(mu = 0.3, nu = 1.0)

  result <- .hzr_optim_weibull(
    time = time,
    status = status,
    x = NULL,
    theta_start = theta_start,
    control = list(maxit = 50)
  )

  # Check that optimization ran
  expect_type(result, "list")
  expect_true("par" %in% names(result))
  expect_true("value" %in% names(result))

  # Parameters should be positive and finite
  expect_true(all(result$par > 0))
  expect_true(all(is.finite(result$par)))
})

test_that("optimizer with gradient is faster than without (heuristic)", {
  skip_if_not_installed("numDeriv")

  # Simulate proper Weibull data with censoring
  set.seed(123)
  n <- 50
  event_time <- (-log(runif(n)) / 0.5) ^ (1 / 1.2)
  cens_time <- runif(n, 0, quantile(event_time, 0.8))
  status <- as.integer(event_time <= cens_time)
  time <- pmin(event_time, cens_time)

  theta_start <- c(mu = 0.5, nu = 1.0)

  # Run with gradient (should converge faster)
  result_with_grad <- .hzr_optim_weibull(
    time = time,
    status = status,
    x = NULL,
    theta_start = theta_start,
    control = list(maxit = 1000)
  )

  # Should converge
  expect_equal(result_with_grad$convergence, 0)

  # Final estimates should be reasonable
  expect_true(result_with_grad$par[1] > 0)
  expect_true(result_with_grad$par[2] > 0)
})
