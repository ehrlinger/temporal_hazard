test_that("hazard() builds a hazard object", {
  x <- matrix(c(1, 0, 0, 1), ncol = 2)
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = x,
    theta = c(0.2, -0.1),
    dist = "weibull",
    maxit = 50
  )

  expect_s3_class(fit, "hazard")
  expect_equal(fit$spec$dist, "weibull")
  expect_equal(fit$legacy_args$maxit, 50)
})

test_that("hazard() validates core dimensions", {
  expect_error(
    hazard(time = c(1, 2), status = 1),
    "same length"
  )

  expect_error(
    hazard(time = c(1, 2), status = c(1, 0), x = matrix(1, nrow = 3, ncol = 1)),
    "rows must match"
  )
})

test_that("predict.hazard returns linear predictor and hazard scale", {
  x <- matrix(c(1, 0, 0, 1), ncol = 2)
  fit <- hazard(time = c(1, 2), status = c(1, 0), x = x, theta = c(0.3, -0.2))

  eta <- predict(fit, type = "linear_predictor")
  hz <- predict(fit, type = "hazard")

  expect_equal(eta, c(0.3, -0.2), tolerance = 1e-12)
  expect_equal(hz, exp(eta), tolerance = 1e-12)
})

test_that("predict.hazard accepts newdata", {
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = matrix(c(1, 2, 3, 4), ncol = 2),
    theta = c(0.5, 0.25)
  )

  newdata <- matrix(c(2, 1, 0, 1), ncol = 2)
  eta <- predict(fit, newdata = newdata, type = "linear_predictor")

  expect_equal(eta, as.numeric(newdata %*% c(0.5, 0.25)), tolerance = 1e-12)
})

test_that("predict.hazard errors when theta is missing", {
  fit <- hazard(time = c(1, 2), status = c(1, 0), x = matrix(c(1, 0), ncol = 1))
  expect_error(predict(fit), "No coefficients")
})

test_that("summary.hazard returns model summary metadata", {
  fit <- hazard(
    time = c(1, 2, 3),
    status = c(1, 0, 1),
    x = matrix(c(1, 0, 0, 1, 1, 1), ncol = 2),
    theta = c(0.3, -0.2),
    dist = "weibull"
  )

  s <- summary(fit)

  expect_s3_class(s, "summary.hazard")
  expect_equal(s$n, 3)
  expect_equal(s$p, 2)
  expect_equal(s$dist, "weibull")
  expect_equal(rownames(s$coefficients), c("beta1", "beta2"))
  expect_equal(s$coefficients$estimate, c(0.3, -0.2), tolerance = 1e-12)
})

test_that("summary.hazard includes standard errors when vcov is available", {
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = matrix(c(1, 0), ncol = 1),
    theta = c(0.4),
    dist = "exponential"
  )
  fit$fit$vcov <- matrix(0.25, nrow = 1, ncol = 1)

  s <- summary(fit)

  expect_true(s$has_vcov)
  expect_equal(s$coefficients$std_error, 0.5, tolerance = 1e-12)
  expect_true(is.finite(s$coefficients$z_stat))
  expect_true(is.finite(s$coefficients$p_value))
})

test_that("print.summary.hazard prints without error", {
  fit <- hazard(
    time = c(1, 2),
    status = c(1, 0),
    x = matrix(c(1, 0), ncol = 1),
    theta = c(0.4),
    dist = "exponential"
  )

  expect_output(print(summary(fit)), "hazard model summary")
})

test_that("hazard() accepts formula interface with Surv()", {
  df <- data.frame(
    time = c(1, 2, 3),
    status = c(1, 0, 1),
    x1 = c(0, 1, 0),
    x2 = c(1, 1, 0)
  )

  fit <- hazard(
    Surv(time, status) ~ x1 + x2,
    data = df,
    theta = c(0.3, -0.2),
    dist = "weibull"
  )

  expect_s3_class(fit, "hazard")
  expect_equal(fit$spec$dist, "weibull")
  expect_equal(length(fit$fit$theta), 2)
})

test_that("formula interface extracts predictors correctly", {
  df <- data.frame(
    time = c(1, 2, 3),
    status = c(1, 0, 1),
    x1 = c(0.5, 1.5, 2.5),
    x2 = c(1, 2, 3)
  )

  fit <- hazard(
    Surv(time, status) ~ x1 + x2,
    data = df,
    theta = c(-0.5, 0.3, -0.1),  # log(lambda), beta1, beta2 for exponential
    dist = "exponential"
  )

  pred <- predict(fit, type = "linear_predictor")
  expected_eta <- c(0.5 * 0.3 + 1 * (-0.1), 1.5 * 0.3 + 2 * (-0.1), 2.5 * 0.3 + 3 * (-0.1))
  expect_equal(pred, expected_eta, tolerance = 1e-12)
})

test_that("formula interface requires data when formula is provided", {
  expect_error(
    hazard(Surv(c(1, 2), c(1, 0)) ~ x, data = NULL),
    "'data' is required"
  )
})

test_that("formula interface works without intercept", {
  df <- data.frame(
    time = c(1, 2, 3),
    status = c(1, 0, 1),
    x = c(1, 0, 1)
  )

  fit <- hazard(
    Surv(time, status) ~ x - 1,
    data = df,
    theta = 0.4,
    dist = "exponential"
  )

  expect_equal(ncol(fit$data$x), 1)
  expect_equal(as.numeric(fit$data$x[, 1]), c(1, 0, 1))
})
