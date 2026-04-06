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
