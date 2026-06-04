# test-weibull-hessian.R -- numDeriv cross-check for .hzr_hessian_weibull_internal()
#
# Reconstruct the optimizer's internal objective (negative log-likelihood on the
# (alpha, psi, beta) scale). Post-PR#58 the natural-scale .hzr_logl_weibull
# equals logl_internal for event+right data, so we can map phi -> (mu, nu, beta)
# and reuse it.
.wb_obj_internal <- function(phi, time, status, time_lower = NULL,
                             time_upper = NULL, x = NULL, weights = NULL) {
  alpha <- phi[1]
  g <- exp(phi[2])
  beta <- if (length(phi) > 2) phi[3:length(phi)] else numeric(0)
  mu <- exp(alpha / g)
  theta_nat <- c(mu, g, beta)
  -.hzr_logl_weibull(theta_nat, time, status, time_lower, time_upper,
                     x = x, weights = weights)
}

test_that(".hzr_hessian_weibull_internal matches numDeriv (no covariates)", {
  skip_if_not_installed("numDeriv")
  set.seed(21)
  n <- 300
  time <- rexp(n, 0.6) + 0.01
  status <- rbinom(n, 1, 0.75)
  phi <- c(alpha = 0.4, psi = log(1.5))
  obj <- function(p) .wb_obj_internal(p, time, status)
  H_an <- .hzr_hessian_weibull_internal(phi, time, status)
  H_nd <- numDeriv::hessian(obj, phi)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_weibull_internal matches numDeriv (covariates + weights)", {
  skip_if_not_installed("numDeriv")
  set.seed(22)
  n <- 400
  x <- cbind(a = rnorm(n), b = rbinom(n, 1, 0.4))
  eta <- as.numeric(x %*% c(0.3, -0.4))
  time <- rweibull(n, shape = 1.4, scale = exp(-eta / 1.4)) + 0.01
  status <- rbinom(n, 1, 0.8)
  w <- runif(n, 0.5, 2)
  phi <- c(alpha = 0.2, psi = log(1.4), 0.3, -0.4)
  obj <- function(p) .wb_obj_internal(p, time, status, x = x, weights = w)
  H_an <- .hzr_hessian_weibull_internal(phi, time, status, x = x, weights = w)
  H_nd <- numDeriv::hessian(obj, phi)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_weibull_internal matches numDeriv (counting-process start)", {
  skip_if_not_installed("numDeriv")
  set.seed(23)
  n <- 250
  start <- runif(n, 0, 0.5)
  time <- start + rexp(n, 0.7) + 0.05
  status <- rbinom(n, 1, 0.7)
  phi <- c(alpha = 0.3, psi = log(1.2))
  obj <- function(p) .wb_obj_internal(p, time, status, time_lower = start)
  H_an <- .hzr_hessian_weibull_internal(phi, time, status, time_lower = start)
  H_nd <- numDeriv::hessian(obj, phi)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_weibull_internal returns NULL for left/interval censoring", {
  set.seed(24)
  n <- 40
  time <- rexp(n, 0.5) + 0.01
  status <- rep(1L, n)
  status[1:5] <- 2L
  expect_null(.hzr_hessian_weibull_internal(
    c(alpha = 0, psi = 0), time, status,
    time_lower = time, time_upper = time + 0.5))
  status2 <- rep(1L, n)
  status2[1:5] <- -1L
  expect_null(.hzr_hessian_weibull_internal(
    c(alpha = 0, psi = 0), time, status2, time_upper = time + 0.5))
})

test_that(".hzr_hessian_weibull_internal handles right-censored time = 0 rows", {
  skip_if_not_installed("numDeriv")
  set.seed(27)
  n <- 60
  time <- c(0, rexp(n - 1, 0.6) + 0.01)      # row 1: right-censored at time 0
  status <- c(0L, rbinom(n - 1, 1, 0.7))
  phi <- c(alpha = 0.3, psi = log(1.3))
  H_an <- .hzr_hessian_weibull_internal(phi, time, status)
  expect_true(all(is.finite(H_an)))
  obj <- function(p) .wb_obj_internal(p, time, status)
  H_nd <- numDeriv::hessian(obj, phi)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})
