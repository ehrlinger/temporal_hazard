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

test_that("summary prints no Hessian note for a clean fit", {
  set.seed(3)
  dat <- data.frame(t = rweibull(250, 1.4, 3), d = rep(1L, 250))
  fit <- hazard(survival::Surv(t, d) ~ 1, data = dat, dist = "weibull",
                fit = TRUE, theta = c(mu = 1, nu = 1))
  out <- capture.output(print(summary(fit)))
  expect_false(any(grepl("ill-conditioned|not positive-definite", out)))
})

test_that("summary prints a note when the fit is flagged ill-conditioned", {
  s <- summary(structure(
    list(spec = list(dist = "weibull", phases = NULL),
         engine = "test",
         fit = list(theta = c(a = 1), vcov = matrix(1, 1, 1),
                    converged = TRUE, objective = -1,
                    rcond = 1e-12, pd = TRUE),
         data = list(time = 1:5, x = NULL),
         call = quote(hazard())),
    class = "hazard"))
  out <- capture.output(print(s))
  expect_true(any(grepl("ill-conditioned", out)))
})

test_that("summary prints a note when the fit is not positive-definite", {
  s <- summary(structure(
    list(spec = list(dist = "weibull", phases = NULL),
         engine = "test",
         fit = list(theta = c(a = 1), vcov = matrix(1, 1, 1),
                    converged = TRUE, objective = -1,
                    rcond = 0.5, pd = FALSE),
         data = list(time = 1:5, x = NULL),
         call = quote(hazard())),
    class = "hazard"))
  out <- capture.output(print(s))
  expect_true(any(grepl("not positive-definite", out)))
  expect_false(any(grepl("ill-conditioned", out)))  # rcond=0.5 is healthy
})
