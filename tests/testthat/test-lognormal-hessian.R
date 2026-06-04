test_that(".hzr_hessian_lognormal matches numDeriv (no covariates)", {
  skip_if_not_installed("numDeriv")
  set.seed(51)
  n <- 300
  time <- exp(rnorm(n, 0.2, 0.8))
  status <- rbinom(n, 1, 0.75)
  theta <- c(mu = 0.2, log_sigma = log(0.8))
  obj <- function(th) -.hzr_logl_lognormal(th, time, status)
  h_an <- .hzr_hessian_lognormal(theta, time, status)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_lognormal matches numDeriv (covariates + weights)", {
  skip_if_not_installed("numDeriv")
  set.seed(52)
  n <- 400
  x <- cbind(a = rnorm(n), b = rbinom(n, 1, 0.4))
  eta <- 0.2 + as.numeric(x %*% c(0.3, -0.4))
  time <- exp(rnorm(n, eta, 0.7))
  status <- rbinom(n, 1, 0.8)
  w <- runif(n, 0.5, 2)
  theta <- c(mu = 0.2, log_sigma = log(0.7), 0.3, -0.4)
  obj <- function(th) -.hzr_logl_lognormal(th, time, status, x = x, weights = w)
  h_an <- .hzr_hessian_lognormal(theta, time, status, x = x, weights = w)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_lognormal matches numDeriv (heavy right-censoring)", {
  skip_if_not_installed("numDeriv")
  set.seed(53)
  n <- 300
  time <- exp(rnorm(n, 0.0, 0.6))
  status <- rbinom(n, 1, 0.3)
  theta <- c(mu = 0.1, log_sigma = log(0.6))
  obj <- function(th) -.hzr_logl_lognormal(th, time, status)
  h_an <- .hzr_hessian_lognormal(theta, time, status)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_lognormal returns NULL for left/interval censoring", {
  set.seed(54)
  n <- 40
  time <- exp(rnorm(n, 0, 0.5))
  status <- rep(1L, n)
  status[1:5] <- 2L
  expect_null(.hzr_hessian_lognormal(
    c(mu = 0, log_sigma = 0), time, status,
    time_lower = time, time_upper = time + 0.5))
  status2 <- rep(1L, n)
  status2[1:5] <- -1L
  expect_null(.hzr_hessian_lognormal(
    c(mu = 0, log_sigma = 0), time, status2, time_upper = time + 0.5))
})

test_that("lognormal fit vcov uses the analytic Hessian (matches numDeriv)", {
  skip_if_not_installed("numDeriv")
  set.seed(55)
  n <- 500
  z <- rnorm(n)
  time <- exp(rnorm(n, 0.2 + 0.5 * z, 0.7))
  status <- rbinom(n, 1, 0.85)
  fit <- hazard(
    survival::Surv(time, status) ~ z,
    data = data.frame(time, status, z),
    dist = "lognormal",
    theta = c(mu = 0, log_sigma = 0, z = 0), fit = TRUE
  )
  obj <- function(th) -.hzr_logl_lognormal(th, time, status, x = cbind(z))
  v_nd <- solve(numDeriv::hessian(obj, fit$fit$theta))
  expect_equal(unname(fit$fit$vcov), unname(v_nd), tolerance = 1e-4)
  expect_true(isTRUE(fit$fit$pd))
})

test_that("lognormal SEs are invariant to covariate rescaling", {
  set.seed(56)
  n <- 600
  z <- rnorm(n)
  time <- exp(rnorm(n, 0.1 + 0.4 * z, 0.6))
  status <- rbinom(n, 1, 0.9)
  f1 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z),
               dist = "lognormal",
               theta = c(mu = 0, log_sigma = 0, z = 0), fit = TRUE)
  f2 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z / 100),
               dist = "lognormal",
               theta = c(mu = 0, log_sigma = 0, z = 0), fit = TRUE)
  se1 <- sqrt(diag(f1$fit$vcov))
  se2 <- sqrt(diag(f2$fit$vcov))
  expect_equal(unname(se1[1]), unname(se2[1]), tolerance = 1e-2)
  expect_equal(unname(se1[2]), unname(se2[2]), tolerance = 1e-2)
})

test_that(".hzr_hessian_lognormal errors on inconsistent theta/x", {
  set.seed(57)
  n <- 30
  time <- exp(rnorm(n, 0, 0.5))
  status <- rbinom(n, 1, 0.7)
  # theta declares one covariate but x is NULL.
  expect_error(
    .hzr_hessian_lognormal(c(mu = 0, log_sigma = 0, z = 0.3), time, status, x = NULL),
    "ncol")
  # x supplied with the wrong number of columns.
  x2 <- cbind(a = rnorm(n), b = rnorm(n))
  expect_error(
    .hzr_hessian_lognormal(c(mu = 0, log_sigma = 0, z = 0.3), time, status, x = x2),
    "ncol")
})
