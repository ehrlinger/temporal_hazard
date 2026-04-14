# tests/testthat/test-diagnostics.R
# Tests for hzr_deciles() and related diagnostic utilities

# Helper cache: fit once to reduce repeated optimizer runs across tests.
.avc_fixture <- local({
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ age + mal,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
    fit   = TRUE
  )
  stopifnot(isTRUE(fit$fit$converged))
  list(avc = avc, fit = fit)
})

.fit_avc_weibull <- function() .avc_fixture$fit
.avc_data <- function() .avc_fixture$avc

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
                      "total_events", "total_expected",
                      "n_included", "n_excluded"))
  expect_equal(ov$time, 120)
  expect_equal(ov$groups, 10)
})

test_that("hzr_deciles group counts sum to evaluable n", {
  fit <- .fit_avc_weibull()
  avc <- .avc_data()
  cal <- hzr_deciles(fit, time = 120)
  included <- with(avc, dead == 1 | int_dead >= 120)
  events_120 <- with(avc, dead == 1 & int_dead <= 120 & included)

  expect_equal(sum(cal$n), sum(included))
  expect_equal(sum(cal$events), sum(events_120))
  expect_equal(attr(cal, "overall")$n_excluded, sum(!included))
})

test_that("hzr_deciles works with non-default groups", {
  fit <- .fit_avc_weibull()
  avc <- .avc_data()
  cal5 <- hzr_deciles(fit, time = 60, groups = 5)
  expect_equal(nrow(cal5), 5)
  expect_equal(sum(cal5$n), sum(avc$dead == 1 | avc$int_dead >= 60))
})

test_that("hzr_deciles events and expected totals increase with horizon", {
  fit <- .fit_avc_weibull()
  avc <- .avc_data()
  status_all <- avc$dead
  time_all <- ifelse(avc$dead == 1, avc$int_dead, 1e9)

  cal60 <- hzr_deciles(fit, time = 60, status = status_all,
                       event_time = time_all)
  cal120 <- hzr_deciles(fit, time = 120, status = status_all,
                        event_time = time_all)

  events_60 <- with(avc, sum(dead == 1 & int_dead <= 60))
  events_120 <- with(avc, sum(dead == 1 & int_dead <= 120))

  expect_equal(sum(cal60$events), events_60)
  expect_equal(sum(cal120$events), events_120)
  expect_equal(sum(cal60$n), nrow(avc))
  expect_equal(sum(cal120$n), nrow(avc))
  expect_true(sum(cal120$events) >= sum(cal60$events))
  expect_true(sum(cal120$expected) >= sum(cal60$expected))
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
  included <- with(cabgkul, dead == 1 | int_dead >= 12)

  expect_s3_class(cal, "hzr_deciles")
  expect_equal(nrow(cal), 10)
  expect_equal(sum(cal$n), sum(included))
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
  included <- with(cabgkul, dead == 1 | int_dead >= 60)

  expect_s3_class(cal, "hzr_deciles")
  expect_equal(nrow(cal), 10)
  expect_equal(sum(cal$n), sum(included))
})

test_that("hzr_deciles rejects invalid inputs", {
  fit <- .fit_avc_weibull()
  avc <- .avc_data()

  expect_error(hzr_deciles(fit, time = -1), "positive")
  expect_error(hzr_deciles(fit, time = 120, groups = 1), "at least 2")
  expect_error(hzr_deciles(fit, time = 120, groups = nrow(avc) + 1),
               "included observations")
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


# =========================================================================
# hzr_gof tests
# =========================================================================

test_that("hzr_gof returns correct structure", {
  fit <- .fit_avc_weibull()
  gof <- hzr_gof(fit)

  expect_s3_class(gof, "hzr_gof")
  expect_s3_class(gof, "data.frame")
  expect_true(nrow(gof) > 0)
  expect_true(all(c("time", "n_risk", "n_event", "n_censor",
                     "km_surv", "km_cumhaz", "par_surv", "par_cumhaz",
                     "cum_observed", "cum_expected", "residual") %in%
                    names(gof)))

  # Summary attribute
  s <- attr(gof, "summary")
  expect_true(is.list(s))
  expect_named(s, c("total_observed", "total_expected", "final_residual",
                     "dist", "n"))
})

test_that("hzr_gof cumulative observed equals total events", {
  fit <- .fit_avc_weibull()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  gof <- hzr_gof(fit)

  s <- attr(gof, "summary")
  expect_equal(s$total_observed, sum(avc$dead))
  expect_equal(s$n, nrow(avc))
})

test_that("hzr_gof survival values are in [0, 1]", {
  fit <- .fit_avc_weibull()
  gof <- hzr_gof(fit)

  expect_true(all(gof$km_surv >= 0 & gof$km_surv <= 1))
  expect_true(all(gof$par_surv >= 0 & gof$par_surv <= 1))
})

test_that("hzr_gof cumulative hazard is non-negative and increasing", {
  fit <- .fit_avc_weibull()
  gof <- hzr_gof(fit)

  expect_true(all(gof$km_cumhaz >= 0))
  expect_true(all(gof$par_cumhaz >= 0))
  # Parametric cumhaz should be non-decreasing
  expect_true(all(diff(gof$par_cumhaz) >= -1e-10))
})

test_that("hzr_gof conservation ratio is reasonable", {
  fit <- .fit_avc_weibull()
  gof <- hzr_gof(fit)

  s <- attr(gof, "summary")
  ratio <- s$total_expected / s$total_observed
  # For a reasonable fit, E/O should be in the ballpark (0.5 to 2.0)
  expect_true(ratio > 0.2 && ratio < 5.0,
              label = paste("conservation ratio", round(ratio, 3),
                            "should be between 0.2 and 5.0"))
})

test_that("hzr_gof works with custom time grid", {
  fit <- .fit_avc_weibull()
  t_grid <- seq(1, 200, by = 10)
  gof <- hzr_gof(fit, time_grid = t_grid)

  expect_equal(nrow(gof), length(t_grid))
  expect_equal(gof$time, t_grid)
})

test_that("hzr_gof works with intercept-only model", {
  data(cabgkul, package = "TemporalHazard")
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = cabgkul,
    dist  = "weibull",
    theta = c(mu = 0.10, nu = 1.0),
    fit   = TRUE
  )
  gof <- hzr_gof(fit)

  expect_s3_class(gof, "hzr_gof")
  expect_true(nrow(gof) > 0)
  s <- attr(gof, "summary")
  expect_equal(s$total_observed, sum(cabgkul$dead))
})

test_that("hzr_gof works with multiphase model", {
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
  gof <- hzr_gof(fit_mp)

  expect_s3_class(gof, "hzr_gof")
  # Should have phase-specific columns
  expect_true(any(grepl("^par_cumhaz_", names(gof))))
})

test_that("hzr_gof rejects invalid inputs", {
  expect_error(hzr_gof("not_a_model"), "hazard object")

  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit_unfitted <- hazard(
    survival::Surv(int_dead, dead) ~ age + mal,
    data  = avc,
    dist  = "weibull",
    theta = c(mu = 0.01, nu = 0.5, 0, 0),
    fit   = FALSE
  )
  expect_error(hzr_gof(fit_unfitted), "fit = TRUE")
})

test_that("print.hzr_gof runs without error", {
  fit <- .fit_avc_weibull()
  gof <- hzr_gof(fit)
  expect_output(print(gof), "Goodness-of-fit")
  expect_output(print(gof), "Conservation ratio")
})


# =========================================================================
# hzr_kaplan tests
# =========================================================================

test_that("hzr_kaplan returns correct structure", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)

  expect_s3_class(km, "hzr_kaplan")
  expect_s3_class(km, "data.frame")
  expect_true(nrow(km) > 0)
  expect_named(km, c("time", "n_risk", "n_event", "n_censor",
                      "survival", "std_err", "cl_lower", "cl_upper",
                      "cumhaz", "hazard", "density", "life"))
})

test_that("hzr_kaplan survival is monotone non-increasing", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)

  expect_true(all(diff(km$survival) <= 0),
              label = "survival should be non-increasing")
  expect_true(all(km$survival >= 0 & km$survival <= 1))
})

test_that("hzr_kaplan confidence limits respect [0, 1]", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)

  expect_true(all(km$cl_lower >= 0 & km$cl_lower <= 1))
  expect_true(all(km$cl_upper >= 0 & km$cl_upper <= 1))
  # Lower <= survival <= upper
  expect_true(all(km$cl_lower <= km$survival + 1e-10))
  expect_true(all(km$cl_upper >= km$survival - 1e-10))
})

test_that("hzr_kaplan matches survfit survival values", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead, event_only = FALSE)
  sf <- survival::survfit(survival::Surv(int_dead, dead) ~ 1,
                           data = cabgkul)

  # Survival values at survfit times should match exactly
  expect_equal(km$survival, sf$surv, tolerance = 1e-12)
  expect_equal(km$time, sf$time)
})

test_that("hzr_kaplan cumhaz is consistent with survival", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)

  # cumhaz = -log(S)
  expected_cumhaz <- -log(pmax(km$survival, .Machine$double.xmin))
  expect_equal(km$cumhaz, expected_cumhaz, tolerance = 1e-10)
})

test_that("hzr_kaplan life integral is non-negative and non-decreasing", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)

  expect_true(all(km$life >= 0))
  expect_true(all(diff(km$life) >= -1e-10),
              label = "life integral should be non-decreasing")
})

test_that("hzr_kaplan event_only filters correctly", {
  data(cabgkul, package = "TemporalHazard")
  km_events <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead,
                            event_only = TRUE)
  km_all <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead,
                         event_only = FALSE)

  expect_true(nrow(km_all) >= nrow(km_events))
  expect_true(all(km_events$n_event > 0))
})

test_that("hzr_kaplan conf_level parameter works", {
  data(cabgkul, package = "TemporalHazard")
  km_95 <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead,
                        conf_level = 0.95)
  km_68 <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead,
                        conf_level = 0.68268948)

  # 95% intervals should be wider than 68% (1-SD) intervals
  width_95 <- km_95$cl_upper - km_95$cl_lower
  width_68 <- km_68$cl_upper - km_68$cl_lower
  # Compare at matching times where both have events
  expect_true(all(width_95 >= width_68 - 1e-10))
})

test_that("hzr_kaplan rejects invalid inputs", {
  expect_error(hzr_kaplan("a", 1), "numeric")
  expect_error(hzr_kaplan(1:5, 1:3), "same length")
  expect_error(hzr_kaplan(1:5, rep(1, 5), conf_level = 1.5),
               "between 0 and 1")
})

test_that("print.hzr_kaplan runs without error", {
  data(cabgkul, package = "TemporalHazard")
  km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)
  expect_output(print(km), "Kaplan-Meier")
  expect_output(print(km), "RMST")
})
