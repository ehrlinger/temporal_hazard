test_that(".hzr_hessian_loglogistic matches numDeriv (no covariates)", {
  skip_if_not_installed("numDeriv")
  set.seed(31)
  n <- 300
  time <- rexp(n, 0.6) + 0.01
  status <- rbinom(n, 1, 0.75)
  theta <- c(log_alpha = log(0.5), log_beta = log(1.4))
  obj <- function(th) -.hzr_logl_loglogistic(th, time, status)
  h_an <- .hzr_hessian_loglogistic(theta, time, status)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_loglogistic matches numDeriv (covariates + weights)", {
  skip_if_not_installed("numDeriv")
  set.seed(32)
  n <- 400
  x <- cbind(a = rnorm(n), b = rbinom(n, 1, 0.4))
  eta <- as.numeric(x %*% c(0.3, -0.4))
  time <- rexp(n, rate = exp(-eta) * 0.5) + 0.01
  status <- rbinom(n, 1, 0.8)
  w <- runif(n, 0.5, 2)
  theta <- c(log_alpha = log(0.5), log_beta = log(1.3), 0.3, -0.4)
  obj <- function(th) -.hzr_logl_loglogistic(th, time, status, x = x, weights = w)
  h_an <- .hzr_hessian_loglogistic(theta, time, status, x = x, weights = w)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_loglogistic handles right-censored time = 0 rows", {
  skip_if_not_installed("numDeriv")
  set.seed(33)
  n <- 60
  time <- c(0, rexp(n - 1, 0.6) + 0.01)
  status <- c(0L, rbinom(n - 1, 1, 0.7))
  theta <- c(log_alpha = log(0.6), log_beta = log(1.2))
  h_an <- .hzr_hessian_loglogistic(theta, time, status)
  expect_true(all(is.finite(h_an)))
  obj <- function(th) -.hzr_logl_loglogistic(th, time, status)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_loglogistic returns NULL for left/interval censoring", {
  set.seed(34)
  n <- 40
  time <- rexp(n, 0.5) + 0.01
  status <- rep(1L, n)
  status[1:5] <- 2L
  expect_null(.hzr_hessian_loglogistic(
    c(log_alpha = 0, log_beta = 0), time, status,
    time_lower = time, time_upper = time + 0.5))
  status2 <- rep(1L, n)
  status2[1:5] <- -1L
  expect_null(.hzr_hessian_loglogistic(
    c(log_alpha = 0, log_beta = 0), time, status2, time_upper = time + 0.5))
})
