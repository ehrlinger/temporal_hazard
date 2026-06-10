# Regression tests for the Weibull event-hazard / cumulative-hazard consistency
# bug: haz_event must equal dH/dt for H = (mu*t)^nu*exp(eta), i.e.
# h = nu*mu^nu*t^(nu-1)*exp(eta) = (nu/t)*H  (Form A, matching C/SAS HAZARD).

test_that(".hzr_logl_weibull event hazard is consistent with its cumulative hazard", {
  set.seed(11)
  n <- 150
  time <- rexp(n, 0.6) + 0.01
  status <- rbinom(n, 1, 0.7)
  mu <- 0.7
  nu <- 1.6
  cumhaz <- (mu * time)^nu                  # Form A cumulative hazard
  haz <- nu * mu^nu * time^(nu - 1)         # correct hazard = dH/dt = (nu/t)*H
  ll_expected <- sum(status * (log(haz) - cumhaz)) - sum((1 - status) * cumhaz)
  ll_actual <- .hzr_logl_weibull(c(mu, nu), time, status)
  expect_equal(ll_actual, ll_expected, tolerance = 1e-8)
})

test_that(".hzr_gradient_weibull matches numDeriv of the (corrected) log-likelihood", {
  skip_if_not_installed("numDeriv")
  set.seed(12)
  n <- 120
  x <- cbind(z = rnorm(n))
  time <- rexp(n, 0.5) + 0.01
  status <- rbinom(n, 1, 0.75)
  theta <- c(mu = 0.8, nu = 1.4, z = 0.3)
  ll <- .hzr_logl_weibull(theta, time, status, x = x, return_gradient = TRUE)
  g_an <- attr(ll, "gradient")
  obj <- function(th) .hzr_logl_weibull(th, time, status, x = x)
  g_nd <- numDeriv::grad(obj, theta)
  expect_equal(g_an, g_nd, tolerance = 1e-5)
})

test_that(".hzr_gradient_weibull is finite for right-censored time = 0 rows", {
  skip_if_not_installed("numDeriv")
  set.seed(41)
  n <- 60
  time <- c(0, rexp(n - 1, 0.5) + 0.05)
  status <- c(0L, rbinom(n - 1, 1, 0.7))
  theta <- c(mu = 0.6, nu = 1.3)
  ll <- .hzr_logl_weibull(theta, time, status, return_gradient = TRUE)
  g_an <- attr(ll, "gradient")
  expect_true(all(is.finite(g_an)))
  obj <- function(th) .hzr_logl_weibull(th, time, status)
  g_nd <- numDeriv::grad(obj, theta)
  expect_equal(g_an, g_nd, tolerance = 1e-5)
})
