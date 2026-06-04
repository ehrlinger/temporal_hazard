test_that(".hzr_hessian_exponential matches numDeriv (no covariates)", {
  skip_if_not_installed("numDeriv")
  set.seed(2)
  n <- 300
  time <- rexp(n, rate = 0.7)
  status <- rbinom(n, 1, 0.7)
  theta <- c(log_rate = log(0.7))
  obj <- function(par) -.hzr_logl_exponential(par, time, status, x = NULL)
  H_an <- .hzr_hessian_exponential(theta, time, status, x = NULL)
  H_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_exponential matches numDeriv (covariates + weights)", {
  skip_if_not_installed("numDeriv")
  set.seed(3)
  n <- 400
  x <- cbind(a = rnorm(n), b = rbinom(n, 1, 0.4))
  eta <- as.numeric(x %*% c(0.3, -0.5))
  time <- rexp(n, rate = exp(log(0.5) + eta))
  status <- rbinom(n, 1, 0.8)
  w <- runif(n, 0.5, 2)
  theta <- c(log_rate = log(0.5), 0.3, -0.5)
  obj <- function(par) -.hzr_logl_exponential(par, time, status, x = x, weights = w)
  H_an <- .hzr_hessian_exponential(theta, time, status, x = x, weights = w)
  H_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_exponential returns NULL for left/interval censoring", {
  set.seed(4)
  n <- 50
  time <- rexp(n, 0.5)
  status <- rep(1L, n)
  status[1:5] <- 2L
  tl <- time
  tu <- time + 0.5
  expect_null(.hzr_hessian_exponential(c(log_rate = 0), time, status,
                                       time_lower = tl, time_upper = tu, x = NULL))
  status2 <- rep(1L, n)
  status2[1:5] <- -1L
  expect_null(.hzr_hessian_exponential(c(log_rate = 0), time, status2,
                                       time_upper = tu, x = NULL))
})

test_that("exponential fit vcov uses the analytic Hessian (matches numDeriv)", {
  skip_if_not_installed("numDeriv")
  set.seed(5)
  n <- 500
  z <- rnorm(n)
  time <- rexp(n, rate = exp(log(0.4) + 0.6 * z))
  status <- rbinom(n, 1, 0.85)
  fit <- hazard(
    survival::Surv(time, status) ~ z,
    data = data.frame(time, status, z),
    dist = "exponential", theta = c(log_rate = 0, z = 0), fit = TRUE
  )
  obj <- function(par) -.hzr_logl_exponential(par, time, status, x = cbind(z))
  v_nd <- solve(numDeriv::hessian(obj, fit$fit$theta))
  expect_equal(unname(fit$fit$vcov), unname(v_nd), tolerance = 1e-4)
  expect_true(isTRUE(fit$fit$pd))
})

test_that("exponential SEs are invariant to covariate rescaling", {
  set.seed(6)
  n <- 600
  z <- rnorm(n, mean = 0, sd = 1)   # centered -> well-conditioned design
  time <- rexp(n, rate = exp(log(0.3) + 0.5 * z))
  status <- rbinom(n, 1, 0.9)
  f1 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z),
               dist = "exponential", theta = c(log_rate = 0, z = 0), fit = TRUE)
  f2 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z / 100),
               dist = "exponential", theta = c(log_rate = 0, z = 0), fit = TRUE)
  se1 <- sqrt(diag(f1$fit$vcov))
  se2 <- sqrt(diag(f2$fit$vcov))
  # log_rate SE is invariant to covariate rescaling.
  expect_equal(unname(se1[1]), unname(se2[1]), tolerance = 1e-3)
  # The beta z-statistic is invariant under x -> x / 100.
  z1 <- unname(f1$fit$theta[2] / se1[2])
  z2 <- unname(f2$fit$theta[2] / se2[2])
  expect_equal(z1, z2, tolerance = 1e-2)
})

test_that(".hzr_hessian_exponential errors on inconsistent theta/x", {
  set.seed(61)
  n <- 30
  time <- rexp(n, 0.5) + 0.01
  status <- rbinom(n, 1, 0.7)
  expect_error(
    .hzr_hessian_exponential(c(log_rate = 0, z = 0.3), time, status, x = NULL),
    "ncol")
  x2 <- cbind(a = rnorm(n), b = rnorm(n))
  expect_error(
    .hzr_hessian_exponential(c(log_rate = 0, z = 0.3), time, status, x = x2),
    "ncol")
})

test_that(".hzr_hessian_exponential guard handles vector x (NCOL robustness)", {
  set.seed(64)
  n <- 30
  time <- rexp(n, 0.5) + 0.01
  status <- rbinom(n, 1, 0.7)
  # x as a bare vector (NCOL = 1) with a 2-covariate theta -> clean "ncol" error,
  # not a cryptic "argument is of length zero" from ncol(vector) == NULL.
  expect_error(
    .hzr_hessian_exponential(c(log_rate = 0, a = 0.1, b = 0.2), time, status,
                             x = rnorm(n)),
    "ncol")
})

test_that(".hzr_hessian_exponential coerces a vector x to a 1-column design", {
  skip_if_not_installed("numDeriv")
  set.seed(65)
  n <- 200
  z <- rnorm(n)
  time <- rexp(n, rate = exp(-0.4 * z) * 0.5) + 0.01
  status <- rbinom(n, 1, 0.8)
  theta <- c(log_rate = log(0.5), z = -0.4)
  h_vec <- .hzr_hessian_exponential(theta, time, status, x = z)            # bare vector
  h_mat <- .hzr_hessian_exponential(theta, time, status, x = matrix(z, n, 1))
  expect_equal(unname(h_vec), unname(h_mat), tolerance = 1e-10)
  expect_equal(dim(h_vec), c(2L, 2L))
})
