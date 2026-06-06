# predict(type = "hazard") for multiphase models: the instantaneous additive
# hazard h(t|x) = sum_j mu_j(x) * phi_j'(t). Single-distribution "hazard"
# (exp(eta)) is unchanged; this covers the multiphase path added so the
# instantaneous hazard has a public predict() route (previously internal only).

.ph_fit <- function(formula = survival::Surv(time, status) ~ 1, data, ...) {
  set.seed(1)
  hazard(formula, data = data, dist = "multiphase",
         phases = list(
           early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                fixed = "shapes", ...),
           constant = hzr_phase("constant")),
         fit = TRUE, control = list(n_starts = 1, maxit = 500, conserve = TRUE))
}

.ph_data <- function(n = 300, seed = 3) {
  set.seed(seed)
  x <- rnorm(n)
  t <- rexp(n, 0.3) + 0.01
  cens <- runif(n, 0.2, 6)
  data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x)
}

test_that("multiphase predict(type = 'hazard') returns the additive hazard", {
  d <- .ph_data()
  fit <- .ph_fit(data = d)
  grid <- data.frame(time = seq(0.05, 5, length.out = 25))

  haz <- predict(fit, newdata = grid, type = "hazard")
  ref <- TemporalHazard:::.hzr_multiphase_hazard(
    grid$time, fit$fit$theta, fit$fit$phases,
    fit$fit$covariate_counts, fit$fit$x_list)

  expect_equal(unname(haz), unname(ref), tolerance = 1e-10)
  expect_true(all(haz >= 0))
  expect_equal(length(haz), nrow(grid))
})

test_that("multiphase hazard requires time and rejects linear_predictor", {
  d <- .ph_data()
  fit <- .ph_fit(data = d)
  expect_error(predict(fit, newdata = data.frame(z = 1), type = "hazard"),
               "time")
  expect_error(predict(fit, type = "linear_predictor"),
               "linear_predictor")
})

test_that("multiphase hazard supports covariate newdata", {
  d <- .ph_data()
  set.seed(1)
  fit <- hazard(survival::Surv(time, status) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1, maxit = 500, conserve = TRUE))
  h0 <- predict(fit, newdata = data.frame(time = 1, x = 0), type = "hazard")
  h1 <- predict(fit, newdata = data.frame(time = 1, x = 1), type = "hazard")
  expect_true(is.finite(h0) && is.finite(h1))
  expect_true(h0 > 0 && h1 > 0)
  expect_false(isTRUE(all.equal(unname(h0), unname(h1))))  # covariate matters
})

test_that("multiphase hazard se.fit returns delta-method limits", {
  skip_if_not_installed("numDeriv")
  d <- .ph_data()
  fit <- .ph_fit(data = d)
  grid <- data.frame(time = c(0.1, 0.5, 1, 3))

  out <- predict(fit, newdata = grid, type = "hazard", se.fit = TRUE)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("fit", "se.fit", "lower", "upper"))
  # Point estimate matches the non-se path.
  pt <- predict(fit, newdata = grid, type = "hazard")
  expect_equal(out$fit, unname(pt), tolerance = 1e-8)
  expect_true(all(is.finite(out$se.fit)) && all(out$se.fit > 0))
  expect_true(all(out$lower <= out$fit + 1e-9) &&
              all(out$upper >= out$fit - 1e-9))
  expect_true(all(out$lower >= 0))  # log-scale CL keeps hazard positive
})

test_that("multiphase hazard rejects decompose = TRUE", {
  d <- .ph_data()
  fit <- .ph_fit(data = d)
  expect_error(
    predict(fit, newdata = data.frame(time = 1), type = "hazard",
            decompose = TRUE),
    "decompose")
})
