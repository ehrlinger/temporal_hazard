# tests/testthat/test-diagnostics.R
# Tests for hzr_deciles() and related diagnostic utilities

# Helper: fit a Weibull model with covariates (produces varying predictions)
.fit_avc_weibull <- function() {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ age + mal,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
    fit   = TRUE
  )
  stopifnot(!is.na(fit$fit$converged))
  fit
}

test_that("hzr_deciles returns correct structure", {
  fit <- .fit_avc_weibull()
  cal <- hzr_deciles(fit, time = 120)

  expect_s3_class(cal, "hzr_deciles")
  expect_s3_class(cal, "data.frame")
  expect_equal(nrow(cal), 10)
  expect_named(cal, c("group", "n", "events", "expected", "observed_rate",
                       "expected_rate", "chi_sq", "p_value",
                       "mean_survival", "mean_cumhaz"))

  # Overall attribute
  ov <- attr(cal, "overall")
  expect_true(is.list(ov))
  expect_named(ov, c("chi_sq", "df", "p_value", "time", "groups",
                      "total_events", "total_expected"))
  expect_equal(ov$time, 120)
  expect_equal(ov$groups, 10)
})

test_that("hzr_deciles group counts sum to n", {
  fit <- .fit_avc_weibull()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  cal <- hzr_deciles(fit, time = 120)

  expect_equal(sum(cal$n), nrow(avc))
  expect_equal(sum(cal$events), sum(avc$dead))
})

test_that("hzr_deciles works with non-default groups", {
  fit <- .fit_avc_weibull()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  cal5 <- hzr_deciles(fit, time = 60, groups = 5)
  expect_equal(nrow(cal5), 5)
  expect_equal(sum(cal5$n), nrow(avc))
})

test_that("hzr_deciles risk ordering is monotone", {
  fit <- .fit_avc_weibull()
  cal <- hzr_deciles(fit, time = 120)

  # Mean cumhaz should increase across groups (group 1 = lowest risk)
  expect_true(all(diff(cal$mean_cumhaz) >= -1e-10),
              label = "mean cumulative hazard should increase across groups")
  # Mean survival should decrease across groups
  expect_true(all(diff(cal$mean_survival) <= 1e-10),
              label = "mean survival should decrease across groups")
})

test_that("hzr_deciles chi-square values are non-negative", {
  fit <- .fit_avc_weibull()
  cal <- hzr_deciles(fit, time = 120)
  valid <- !is.na(cal$chi_sq)
  expect_true(all(cal$chi_sq[valid] >= 0))
  expect_true(all(cal$p_value[valid] >= 0 & cal$p_value[valid] <= 1))

  ov <- attr(cal, "overall")
  expect_true(ov$chi_sq >= 0)
  expect_true(ov$df > 0)
})

test_that("hzr_deciles works with intercept-only model", {
  # Intercept-only: all predictions identical, but grouping should still work
  data(cabgkul, package = "TemporalHazard")
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = cabgkul,
    dist  = "weibull",
    theta = c(mu = 0.10, nu = 1.0),
    fit   = TRUE
  )
  cal <- hzr_deciles(fit, time = 12)

  expect_s3_class(cal, "hzr_deciles")
  expect_equal(nrow(cal), 10)
  expect_equal(sum(cal$n), nrow(cabgkul))
  # All expected rates should be equal (same prediction for everyone)
  expect_true(max(cal$expected_rate) - min(cal$expected_rate) < 1e-10)
})

test_that("hzr_deciles works with multiphase model", {
  data(cabgkul, package = "TemporalHazard")
  fit_mp <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = cabgkul,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1,
                            fixed = "shapes"),
      constant = hzr_phase("constant"),
      late     = hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1,
                            fixed = "shapes")
    ),
    fit     = TRUE,
    control = list(n_starts = 3, maxit = 500)
  )
  cal <- hzr_deciles(fit_mp, time = 60)

  expect_s3_class(cal, "hzr_deciles")
  expect_equal(nrow(cal), 10)
  expect_equal(sum(cal$n), nrow(cabgkul))
})

test_that("hzr_deciles rejects invalid inputs", {
  fit <- .fit_avc_weibull()

  expect_error(hzr_deciles(fit, time = -1), "positive")
  expect_error(hzr_deciles(fit, time = 120, groups = 1), "at least 2")
  expect_error(hzr_deciles("not_a_model", time = 12), "hazard object")
})

test_that("hzr_deciles unfitted model is rejected", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit_unfitted <- hazard(
    survival::Surv(int_dead, dead) ~ age + status,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.05, nu = 0.5, 0, 0),
    fit   = FALSE
  )

  expect_error(hzr_deciles(fit_unfitted, time = 120), "fit = TRUE")
})

test_that("print.hzr_deciles runs without error", {
  fit <- .fit_avc_weibull()
  cal <- hzr_deciles(fit, time = 120)
  expect_output(print(cal), "Decile-of-risk calibration")
  expect_output(print(cal), "Overall")
})
