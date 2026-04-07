library(testthat)
library(devtools)

load_all()

.make_mixed_censor_data <- function(n = 70, seed = 321) {
  set.seed(seed)
  base_time <- rexp(n, rate = 0.4) + 0.05
  status <- sample(c(1, 0, -1, 2), n, replace = TRUE, prob = c(0.5, 0.3, 0.1, 0.1))

  time <- base_time
  time_lower <- base_time
  time_upper <- base_time

  idx_left <- which(status == -1)
  if (length(idx_left) > 0) {
    time_upper[idx_left] <- base_time[idx_left] + 0.2
    time_lower[idx_left] <- pmax(base_time[idx_left] - 0.2, 0.01)
    time[idx_left] <- time_upper[idx_left]
  }

  idx_interval <- which(status == 2)
  if (length(idx_interval) > 0) {
    time_lower[idx_interval] <- pmax(base_time[idx_interval] * 0.8, 0.01)
    time_upper[idx_interval] <- base_time[idx_interval] * 1.3 + 0.05
    time[idx_interval] <- (time_lower[idx_interval] + time_upper[idx_interval]) / 2
  }

  list(time = time, status = status, time_lower = time_lower, time_upper = time_upper)
}

test_that("Exponential log-likelihood supports mixed censoring", {
  d <- .make_mixed_censor_data(50, 11)
  ll <- .hzr_logl_exponential(
    theta = c(log_lambda = log(0.4)),
    time = d$time,
    status = d$status,
    time_lower = d$time_lower,
    time_upper = d$time_upper
  )
  expect_true(is.finite(ll))
})

test_that("Log-logistic log-likelihood supports mixed censoring", {
  d <- .make_mixed_censor_data(50, 12)
  ll <- .hzr_logl_loglogistic(
    theta = c(log_alpha = log(0.6), log_beta = log(1.3)),
    time = d$time,
    status = d$status,
    time_lower = d$time_lower,
    time_upper = d$time_upper
  )
  expect_true(is.finite(ll))
})

test_that("Log-normal log-likelihood supports mixed censoring", {
  d <- .make_mixed_censor_data(50, 13)
  ll <- .hzr_logl_lognormal(
    theta = c(mu = 0.8, log_sigma = log(0.7)),
    time = d$time,
    status = d$status,
    time_lower = d$time_lower,
    time_upper = d$time_upper
  )
  expect_true(is.finite(ll))
})

test_that("hazard() fits exponential with mixed censoring", {
  d <- .make_mixed_censor_data(70, 21)
  fit <- hazard(
    time = d$time,
    status = d$status,
    time_lower = d$time_lower,
    time_upper = d$time_upper,
    theta = c(log_lambda = log(0.5)),
    dist = "exponential",
    fit = TRUE,
    control = list(maxit = 250)
  )

  expect_true(is.finite(fit$fit$objective))
  expect_equal(length(fit$fit$theta), 1)
})

test_that("hazard() fits log-logistic with mixed censoring", {
  d <- .make_mixed_censor_data(70, 22)
  fit <- hazard(
    time = d$time,
    status = d$status,
    time_lower = d$time_lower,
    time_upper = d$time_upper,
    theta = c(log_alpha = log(0.6), log_beta = log(1.2)),
    dist = "loglogistic",
    fit = TRUE,
    control = list(maxit = 300)
  )

  expect_true(is.finite(fit$fit$objective))
  expect_equal(length(fit$fit$theta), 2)
})

test_that("hazard() fits log-normal with mixed censoring", {
  d <- .make_mixed_censor_data(70, 23)
  fit <- hazard(
    time = d$time,
    status = d$status,
    time_lower = d$time_lower,
    time_upper = d$time_upper,
    theta = c(mu = 0.8, log_sigma = log(0.9)),
    dist = "lognormal",
    fit = TRUE,
    control = list(maxit = 300)
  )

  expect_true(is.finite(fit$fit$objective))
  expect_equal(length(fit$fit$theta), 2)
})

test_that("Interval bounds validation is enforced for all distributions", {
  bad_lower <- c(1.5, 2.0)
  bad_upper <- c(1.0, 2.5)
  time <- c(1.2, 2.2)
  status <- c(2, 1)

  expect_true(is.infinite(.hzr_logl_exponential(
    theta = c(log(0.3)),
    time = time,
    status = status,
    time_lower = bad_lower,
    time_upper = bad_upper
  )))

  expect_true(is.infinite(.hzr_logl_loglogistic(
    theta = c(log(0.6), log(1.2)),
    time = time,
    status = status,
    time_lower = bad_lower,
    time_upper = bad_upper
  )))

  expect_true(is.infinite(.hzr_logl_lognormal(
    theta = c(0.8, log(0.7)),
    time = time,
    status = status,
    time_lower = bad_lower,
    time_upper = bad_upper
  )))
})
