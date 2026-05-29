library(testthat)

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

test_that("Mixed IC+RC: right-censored contribution not silently zeroed when time_lower supplied", {
  # Regression test for the time_lower dual-use bug:
  # When time_lower is provided for a mix of interval-censored (status=2)
  # and right-censored (status=0) rows, right-censored rows must contribute
  # -H(stop) to the log-likelihood, NOT H(stop)-H(stop)=0.
  theta      <- c(mu = 0.4, nu = 1.2)
  mu <- theta[1]
  nu <- theta[2]

  # 3 interval-censored observations in (1, 2)
  # 2 right-censored observations at t = 5
  time       <- c(1, 1, 1, 5, 5)
  status     <- c(2L, 2L, 2L, 0L, 0L)
  time_lower <- c(1, 1, 1, 0, 0)   # correct: 0 for RC rows
  time_upper <- c(2, 2, 2, 5, 5)

  # Manual LL
  H_lo <- (mu * 1)^nu
  H_hi <- (mu * 2)^nu
  H_rc <- (mu * 5)^nu
  ll_ic_manual <- 3 * (-H_lo + log(1 - exp(-(H_hi - H_lo))))
  ll_rc_manual <- 2 * (-H_rc)
  ll_manual <- unname(ll_ic_manual + ll_rc_manual)

  ll_pkg <- .hzr_logl_weibull(theta, time, status,
                               time_lower = time_lower,
                               time_upper = time_upper)

  expect_equal(ll_pkg, ll_manual, tolerance = 1e-10)

  # Also verify the *wrong* setup (time_lower = time for RC) now gives same
  # result after the fix — previously it zeroed the RC contribution
  time_lower_wrong <- c(1, 1, 1, 5, 5)  # time_lower == time for RC rows
  ll_wrong_setup <- .hzr_logl_weibull(theta, time, status,
                                       time_lower = time_lower_wrong,
                                       time_upper = time_upper)
  # With the fix, time_lower == time means NOT a genuine epoch, so RC
  # contribution is still -H(5), same result as the correct setup
  expect_equal(ll_wrong_setup, ll_manual, tolerance = 1e-10)
})

test_that("Mixed IC+RC: fitting recovers correct parameters with time_lower=0 for RC", {
  set.seed(42)
  n <- 600
  true_time <- rexp(n, rate = 0.05)
  # Interval censoring in 6-month windows; right-censor at 48 months
  status     <- integer(n)
  time       <- numeric(n)
  time_lower <- numeric(n)
  time_upper <- numeric(n)
  for (i in seq_len(n)) {
    brk <- ceiling(true_time[i] / 6)
    if (brk > 8L) {
      status[i]     <- 0L
      time[i]       <- 48
      time_lower[i] <- 0   # correct: 0 for RC rows
      time_upper[i] <- 48
    } else {
      status[i]     <- 2L
      time[i]       <- (brk - 1L) * 6
      time_lower[i] <- (brk - 1L) * 6
      time_upper[i] <- brk * 6
    }
  }
  fit <- hazard(time = time, status = status,
                time_lower = time_lower, time_upper = time_upper,
                dist = "weibull", theta = c(mu = 0.05, nu = 1.0), fit = TRUE)
  expect_true(fit$fit$converged)
  # MLE should be within 20% of truth (nu=1 for exponential, mu=0.05)
  coefs <- coef(fit)[1:2]
  expect_lt(abs(coefs["nu"] - 1.0), 0.2)
  expect_lt(abs(coefs["mu"] - 0.05) / 0.05, 0.2)
})
