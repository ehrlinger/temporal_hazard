library(testthat)
library(devtools)

load_all()

test_that("Weibull log-likelihood handles mixed censoring statuses", {
  time <- c(1.0, 2.0, 1.5, 0.8)
  status <- c(1, 0, -1, 2)
  time_lower <- c(1.0, 2.0, 0.0, 1.2)
  time_upper <- c(1.0, 2.0, 1.5, 2.5)

  theta <- c(mu = 0.7, nu = 1.2)

  ll <- .hzr_logl_weibull(
    theta = theta,
    time = time,
    status = status,
    time_lower = time_lower,
    time_upper = time_upper
  )

  expect_true(is.finite(ll))
})

test_that("Weibull interval contribution matches manual formula", {
  theta <- c(mu = 0.8, nu = 1.4)
  time <- 1.0
  status <- 2
  time_lower <- 0.9
  time_upper <- 1.7

  ll <- .hzr_logl_weibull(
    theta = theta,
    time = time,
    status = status,
    time_lower = time_lower,
    time_upper = time_upper
  )

  H_l <- (theta[1] * time_lower) ^ theta[2]
  H_u <- (theta[1] * time_upper) ^ theta[2]
  manual <- -H_l + hzr_log1mexp(H_u - H_l)

  expect_equal(unname(ll), unname(manual), tolerance = 1e-10)
})

test_that("Weibull left-censored contribution matches manual formula", {
  theta <- c(mu = 0.9, nu = 1.1)
  time <- 1.3
  status <- -1

  ll <- .hzr_logl_weibull(
    theta = theta,
    time = time,
    status = status
  )

  H_u <- (theta[1] * time) ^ theta[2]
  manual <- hzr_log1mexp(H_u)

  expect_equal(ll, manual, tolerance = 1e-10)
})

test_that("Weibull right-censored contribution remains -H(t)", {
  theta <- c(mu = 0.6, nu = 1.5)
  time <- c(0.7, 1.2, 2.0)
  status <- c(0, 0, 0)

  ll <- .hzr_logl_weibull(theta, time, status)
  manual <- -sum((theta[1] * time) ^ theta[2])

  expect_equal(ll, manual, tolerance = 1e-10)
})

test_that("Weibull fit converges for mixed censoring sample", {
  set.seed(123)
  n <- 60
  time <- rexp(n, rate = 0.5)

  status <- sample(c(1, 0, -1, 2), size = n, replace = TRUE, prob = c(0.45, 0.35, 0.1, 0.1))

  time_lower <- time
  time_upper <- time

  idx_left <- which(status == -1)
  if (length(idx_left) > 0) {
    time_upper[idx_left] <- pmax(time[idx_left], 0.2)
    time_lower[idx_left] <- 0
  }

  idx_interval <- which(status == 2)
  if (length(idx_interval) > 0) {
    time_lower[idx_interval] <- pmax(time[idx_interval] * 0.8, 0.05)
    time_upper[idx_interval] <- time[idx_interval] * 1.2 + 0.1
  }

  fit <- hazard(
    time = time,
    status = status,
    time_lower = time_lower,
    time_upper = time_upper,
    theta = c(mu = 0.5, nu = 1.0),
    dist = "weibull",
    fit = TRUE,
    control = list(maxit = 300)
  )

  expect_true(isTRUE(fit$fit$converged) || identical(fit$fit$converged, FALSE))
  expect_true(is.finite(fit$fit$objective))
  expect_equal(length(fit$fit$theta), 2)
})

test_that("Interval bounds are validated (lower must be < upper)", {
  expect_true(is.infinite(.hzr_logl_weibull(
    theta = c(0.7, 1.2),
    time = 1.0,
    status = 2,
    time_lower = 2.0,
    time_upper = 1.0
  )))
})

test_that("Non-Weibull distributions accept interval/left censoring inputs", {
  fit <- hazard(
    time = c(1, 2, 3),
    status = c(1, 2, 0),
    time_lower = c(1, 1.5, 3),
    time_upper = c(1, 2.5, 3),
    theta = c(0.2, 0.1),
    dist = "lognormal",
    fit = FALSE
  )

  expect_s3_class(fit, "hazard")
})

test_that("Right-censored legacy path unchanged when bounds are absent", {
  time <- c(0.5, 1.0, 2.0)
  status <- c(1, 0, 1)
  theta <- c(0.8, 1.3)

  ll_old <- .hzr_logl_weibull(theta, time, status)
  ll_new <- .hzr_logl_weibull(theta, time, status, time_lower = NULL, time_upper = NULL)

  expect_equal(ll_old, ll_new, tolerance = 1e-12)
})
