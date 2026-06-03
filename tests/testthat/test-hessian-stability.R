test_that("fit object carries rcond and pd diagnostics (weibull)", {
  set.seed(2)
  n <- 250L
  dat <- data.frame(t = rweibull(n, shape = 1.4, scale = 3),
                    d = rep(1L, n))
  fit <- hazard(survival::Surv(t, d) ~ 1, data = dat, dist = "weibull",
                theta = c(mu = 1, nu = 1), fit = TRUE)
  expect_true(is.finite(fit$fit$rcond))
  expect_true(isTRUE(fit$fit$pd))
})

test_that("fit object carries rcond and pd diagnostics (multiphase)", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                    fixed = "m"),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 3, maxit = 500, conserve = TRUE)
  )
  expect_true(is.finite(fit$fit$rcond))
  expect_false(is.null(fit$fit$pd))
})
