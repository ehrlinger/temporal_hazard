# tests/testthat/test-predict-cl.R
# Phase 4g -- predict.hazard(..., se.fit = TRUE, level = 0.95)
#
# Coverage:
# (1) API contract: output shape, level validation, decompose incompat
# (2) Jacobian parity: analytic vs numeric for Weibull + multiphase
# (3) Delta-method scale: log-CL for hazard/cumhaz, log-log for survival
# (4) Monotonicity: lower <= fit <= upper
# (5) Survival CLs stay in [0, 1]; cumhaz CLs stay in [0, Inf)
# (6) `vcov` NA (non-converged fit): returns NA CLs with a warning
# (7) Fixed shapes / CoE: free-parameter submatrix recovers meaningful CLs
# (8) Default of `se.fit = FALSE` is byte-identical to pre-0.9.8 predict()

make_toy <- function(seed = 17, n = 60) {
  set.seed(seed)
  data.frame(
    time = rexp(n, 0.3),
    status = rbinom(n, 1, 0.6),
    x = rnorm(n)
  )
}

# ---------------------------------------------------------------------------
# (1) API contract
# ---------------------------------------------------------------------------

test_that("se.fit returns a data frame with fit/se.fit/lower/upper", {
  skip_on_cran()
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)
  res <- predict(fit, newdata = data.frame(time = c(0.5, 1, 2), x = 0),
                 type = "survival", se.fit = TRUE)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("fit", "se.fit", "lower", "upper"))
  expect_equal(nrow(res), 3L)
  expect_true(all(is.finite(res$fit)))
  expect_true(all(is.finite(res$se.fit)))
  expect_true(all(res$se.fit >= 0))
})

test_that("invalid level triggers a clean error", {
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)
  expect_error(
    predict(fit, newdata = data.frame(time = 1, x = 0),
            type = "survival", se.fit = TRUE, level = 1.1),
    "level"
  )
  expect_error(
    predict(fit, newdata = data.frame(time = 1, x = 0),
            type = "survival", se.fit = TRUE, level = 0),
    "level"
  )
})

test_that("se.fit = TRUE with decompose = TRUE errors clearly", {
  df <- make_toy()
  phases <- list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                     fixed = "shapes"),
                 constant = hzr_phase("constant"))
  fit <- hazard(survival::Surv(time, status) ~ 1, data = df,
                dist = "multiphase", phases = phases, fit = TRUE,
                control = list(n_starts = 2))
  expect_error(
    predict(fit, newdata = data.frame(time = c(0.5, 1)),
            type = "cumulative_hazard", se.fit = TRUE, decompose = TRUE),
    "decompose"
  )
})

# ---------------------------------------------------------------------------
# (2) Jacobian parity: analytic vs numeric
# ---------------------------------------------------------------------------

test_that("Weibull analytic jacobian matches numDeriv (cumhaz)", {
  skip_if_not_installed("numDeriv")
  set.seed(3)
  t_new <- c(0.3, 0.8, 1.5, 2.2)
  x_new <- matrix(c(-0.5, 0, 0.4, 1.2), ncol = 1)
  theta <- c(mu = 0.6, nu = 1.2, beta = 0.25)

  cumhaz_fn <- function(th) {
    (th[1] * t_new) ^ th[2] * exp(as.numeric(x_new %*% th[3:length(th)]))
  }
  J_ana <- TemporalHazard:::.hzr_predict_jacobian_weibull(
    "cumulative_hazard", theta, t_new, x_new, length(theta)
  )
  J_num <- numDeriv::jacobian(cumhaz_fn, theta)
  expect_equal(J_ana, J_num, tolerance = 1e-6)
})

test_that("Weibull analytic jacobian matches numDeriv (hazard / relative)", {
  skip_if_not_installed("numDeriv")
  x_new <- matrix(c(-0.5, 0, 0.4), ncol = 1)
  theta <- c(mu = 0.6, nu = 1.2, beta = 0.25)
  haz_fn <- function(th) exp(as.numeric(x_new %*% th[3:length(th)]))
  J_ana <- TemporalHazard:::.hzr_predict_jacobian_weibull(
    "hazard", theta, time = NULL, x_new, length(theta)
  )
  J_num <- numDeriv::jacobian(haz_fn, theta)
  expect_equal(J_ana, J_num, tolerance = 1e-6)
})

test_that("Multiphase analytic jacobian matches numDeriv (cumhaz)", {
  skip_if_not_installed("numDeriv")
  skip_on_cran()
  t_new <- c(0.3, 0.8, 1.5, 2.2)
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  cov_counts <- c(early = 0L, constant = 0L)
  x_list <- list(early = NULL, constant = NULL)
  theta <- c(-3, log(0.3), 1, 1, -2)

  J_ana <- TemporalHazard:::.hzr_predict_jacobian_multiphase(
    theta, t_new, phases, cov_counts, x_list, length(theta)
  )
  J_num <- numDeriv::jacobian(
    function(th) {
      TemporalHazard:::.hzr_multiphase_cumhaz(
        t_new, th, phases, cov_counts, x_list
      )
    },
    theta
  )
  expect_equal(J_ana, J_num, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# (3, 4, 5) Scale, monotonicity, and range
# ---------------------------------------------------------------------------

test_that("survival CLs stay in [0, 1] and straddle the point estimate", {
  skip_on_cran()
  df <- make_toy()
  for (dist in c("weibull", "exponential", "loglogistic", "lognormal")) {
    theta_init <- switch(
      dist,
      weibull = c(0.5, 1, 0),
      exponential = c(log(0.3), 0),
      loglogistic = c(0, 0, 0),
      lognormal = c(0, 0, 0)
    )
    fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = dist,
                  theta = theta_init, fit = TRUE)
    res <- predict(fit, newdata = data.frame(time = c(0.5, 1, 2), x = 0),
                   type = "survival", se.fit = TRUE)
    expect_true(all(res$lower >= 0), info = dist)
    expect_true(all(res$upper <= 1), info = dist)
    expect_true(all(res$lower <= res$fit + 1e-8), info = dist)
    expect_true(all(res$fit <= res$upper + 1e-8), info = dist)
  }
})

test_that("cumulative_hazard CLs are strictly positive and log-symmetric", {
  skip_on_cran()
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)
  res <- predict(fit, newdata = data.frame(time = c(0.5, 1, 2), x = 0),
                 type = "cumulative_hazard", se.fit = TRUE)
  expect_true(all(res$lower > 0))
  expect_true(all(res$upper > res$lower))
  # log-symmetric: log(upper) - log(fit) ~= log(fit) - log(lower)
  log_halfwidth_upper <- log(res$upper) - log(res$fit)
  log_halfwidth_lower <- log(res$fit) - log(res$lower)
  expect_equal(log_halfwidth_upper, log_halfwidth_lower, tolerance = 1e-6)
})

test_that("linear_predictor CLs are symmetric on the natural scale", {
  skip_on_cran()
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)
  res <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)),
                 type = "linear_predictor", se.fit = TRUE)
  halfwidth_upper <- res$upper - res$fit
  halfwidth_lower <- res$fit - res$lower
  expect_equal(halfwidth_upper, halfwidth_lower, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# (6) Missing vcov -> NA CLs with a warning
# ---------------------------------------------------------------------------

test_that("missing vcov produces NA CLs with a warning", {
  skip_on_cran()
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)
  # Knock out vcov to simulate a non-converged / ill-conditioned fit
  fit$fit$vcov <- NULL
  expect_warning(
    res <- predict(fit, newdata = data.frame(time = 1, x = 0),
                   type = "survival", se.fit = TRUE),
    "Variance-covariance"
  )
  expect_true(all(is.na(res$se.fit)))
  expect_true(all(is.na(res$lower)))
  expect_true(all(is.na(res$upper)))
  expect_false(anyNA(res$fit))
})

# ---------------------------------------------------------------------------
# (7) Fixed-shapes multiphase recovers CLs from the free-parameter submatrix
# ---------------------------------------------------------------------------

test_that("multiphase CLs work with fixed shapes + CoE (free submatrix)", {
  skip_on_cran()
  df <- make_toy()
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(survival::Surv(time, status) ~ 1, data = df,
                dist = "multiphase", phases = phases, fit = TRUE,
                control = list(n_starts = 2))
  res <- predict(fit, newdata = data.frame(time = c(0.5, 1, 2)),
                 type = "survival", se.fit = TRUE)
  expect_true(all(is.finite(res$se.fit)))
  expect_true(all(res$se.fit > 0))
  expect_true(all(res$lower >= 0))
  expect_true(all(res$upper <= 1))
})

# ---------------------------------------------------------------------------
# (8) Backward compat: se.fit = FALSE reproduces the old scalar-vector return
# ---------------------------------------------------------------------------

test_that("se.fit = FALSE preserves pre-0.9.8 return shape", {
  skip_on_cran()
  df <- make_toy()
  fit <- hazard(survival::Surv(time, status) ~ x, data = df, dist = "weibull",
                theta = c(0.5, 1, 0), fit = TRUE)

  # Each type should return a bare numeric vector with se.fit = FALSE.
  for (type in c("survival", "cumulative_hazard", "hazard",
                 "linear_predictor")) {
    nd <- if (type %in% c("survival", "cumulative_hazard")) {
      data.frame(time = c(0.5, 1, 2), x = 0)
    } else {
      data.frame(x = c(-1, 0, 1))
    }
    out <- predict(fit, newdata = nd, type = type, se.fit = FALSE)
    expect_type(out, "double")
    expect_length(out, 3L)
    expect_false(is.data.frame(out))
  }
})
