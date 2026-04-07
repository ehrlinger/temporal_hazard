library(testthat)

test_that("log-logistic fixture refits reproducibly", {
  fixture_file <- system.file("fixtures", "hz_loglogistic.rds", package = "TemporalHazard")
  skip_if_not(file.exists(fixture_file))

  fixture <- readRDS(fixture_file)

  refit <- hazard(
    time = fixture$data$time,
    status = fixture$data$status,
    theta = c(log_alpha = 0.0, log_beta = 0.2),
    dist = "loglogistic",
    fit = TRUE,
    control = list(maxit = 500, reltol = 1e-6)
  )

  expect_true(is.finite(refit$fit$objective))
  expect_true(all(is.finite(refit$fit$theta)))
  expect_equal(refit$fit$objective, fixture$fit$logl, tolerance = 5e-2)
})

test_that("log-normal fixture refits reproducibly", {
  fixture_file <- system.file("fixtures", "hz_lognormal.rds", package = "TemporalHazard")
  skip_if_not(file.exists(fixture_file))

  fixture <- readRDS(fixture_file)

  refit <- hazard(
    time = fixture$data$time,
    status = fixture$data$status,
    theta = c(mu = 1.0, log_sigma = 0.0),
    dist = "lognormal",
    fit = TRUE,
    control = list(maxit = 500, reltol = 1e-6)
  )

  expect_true(is.finite(refit$fit$objective))
  expect_true(all(is.finite(refit$fit$theta)))
  expect_equal(refit$fit$objective, fixture$fit$logl, tolerance = 1e-6)
})

test_that("near-zero positive times remain numerically stable across distributions", {
  time <- c(1e-6, 2e-6, 5e-6, 1e-5, 2e-5)
  status <- c(1, 0, 1, 0, 1)

  expect_true(is.finite(.hzr_logl_weibull(c(0.7, 1.2), time, status)))
  expect_true(is.finite(.hzr_logl_exponential(c(log(0.3)), time, status)))
  expect_true(is.finite(.hzr_logl_loglogistic(c(log(0.7), log(1.2)), time, status)))
  expect_true(is.finite(.hzr_logl_lognormal(c(0.8, log(0.6)), time, status)))
})

test_that("high-dimensional (p > n) fit returns finite objective", {
  set.seed(2026)
  n <- 25
  p <- 40

  time <- rexp(n, rate = 0.5) + 0.05
  status <- rbinom(n, 1, 0.6)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)

  theta <- c(log(0.4), rep(0, p))

  fit <- hazard(
    time = time,
    status = status,
    x = x,
    theta = theta,
    dist = "exponential",
    fit = TRUE,
    control = list(maxit = 120)
  )

  expect_true(is.finite(fit$fit$objective))
  expect_equal(length(fit$fit$theta), p + 1)
})
