# Shared test fixture: formula-based single-distribution fit with a
# known significant covariate and a noise covariate in `data`.
.fit_stepwise_base <- function(n = 400L, seed = 101L,
                                signal_beta = 0.8) {
  set.seed(seed)
  df <- data.frame(
    time   = rexp(n, rate = 1),
    status = rep(1L, n),
    x1     = rnorm(n),   # candidate in data, strong effect
    x2     = rnorm(n),   # candidate in data, no effect
    x3     = rnorm(n)    # candidate in data, no effect
  )
  # Apply signal to x1
  df$time <- df$time * exp(-signal_beta * df$x1)

  fit <- hazard(
    Surv(time, status) ~ 1,
    data = df,
    theta = c(0.5, 1.0),
    dist = "weibull",
    fit = TRUE
  )
  list(fit = fit, data = df)
}


# .hzr_refit_with_scope -----------------------------------------------------

test_that("refit adds a covariate and grows theta", {
  obj <- .fit_stepwise_base()
  expect_length(obj$fit$fit$theta, 2L)

  refit <- .hzr_refit_with_scope(
    obj$fit, action = "add", var = "x1",
    data = obj$data
  )
  expect_s3_class(refit, "hazard")
  expect_length(refit$fit$theta, 3L)
  expect_true(refit$fit$converged)
  # The signal covariate should have a non-zero estimate
  last_coef <- refit$fit$theta[length(refit$fit$theta)]
  expect_gt(abs(last_coef), 0.1)
})

test_that("refit drops a covariate and shrinks theta", {
  obj <- .fit_stepwise_base()
  with_x <- .hzr_refit_with_scope(obj$fit, action = "add", var = "x1",
                                    data = obj$data)
  back <- .hzr_refit_with_scope(with_x, action = "drop", var = "x1",
                                 data = obj$data)
  expect_length(back$fit$theta, 2L)
  expect_true(back$fit$converged)
})

test_that("refit errors without a formula-based base fit", {
  fit <- .fit_weibull(n = 100L, seed = 1L)
  expect_error(
    .hzr_refit_with_scope(fit, action = "add", var = "x1",
                           data = data.frame()),
    "formula interface"
  )
})

test_that("refit rejects missing `phase` on multiphase fits", {
  data(avc)
  set.seed(1)
  fit <- hazard(
    Surv(int_dead, dead) ~ 1, data = avc,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE,
    control = list(n_starts = 2L, maxit = 500L)
  )
  expect_error(
    .hzr_refit_with_scope(fit, action = "add", var = "age",
                           data = avc),
    "`phase` is required"
  )
})

test_that("refit adds a var to a specific phase for multiphase", {
  data(avc)
  set.seed(1)
  fit <- hazard(
    Surv(int_dead, dead) ~ 1, data = avc,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE,
    control = list(n_starts = 2L, maxit = 500L)
  )
  refit <- .hzr_refit_with_scope(
    fit, action = "add", var = "age", phase = "early",
    data = avc, control = list(n_starts = 2L, maxit = 500L)
  )
  expect_s3_class(refit, "hazard")
  expect_true("early.age" %in% names(coef(refit)))
  expect_false("constant.age" %in% names(coef(refit)))
})


# .hzr_stepwise_candidates --------------------------------------------------

test_that("candidates excludes already-in-model vars and force_out", {
  obj <- .fit_stepwise_base()
  c1 <- .hzr_stepwise_candidates(obj$fit, scope = NULL,
                                   data = obj$data)
  vars <- vapply(c1, `[[`, character(1L), "var")
  expect_setequal(vars, c("x1", "x2", "x3"))

  with_x1 <- .hzr_refit_with_scope(obj$fit, action = "add", var = "x1",
                                     data = obj$data)
  c2 <- .hzr_stepwise_candidates(with_x1, scope = NULL,
                                   data = obj$data,
                                   force_out = "x3")
  vars2 <- vapply(c2, `[[`, character(1L), "var")
  expect_setequal(vars2, "x2")
})

test_that("candidates for multiphase emits (var, phase) pairs", {
  data(avc)
  set.seed(1)
  fit <- hazard(
    Surv(int_dead, dead) ~ 1, data = avc,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE,
    control = list(n_starts = 2L, maxit = 500L)
  )
  scope <- list(
    early    = ~ age + mal,
    constant = ~ age
  )
  cands <- .hzr_stepwise_candidates(fit, scope = scope, data = avc)
  pairs <- vapply(cands, function(c) paste0(c$var, "@", c$phase),
                  character(1L))
  expect_setequal(pairs, c("age@early", "mal@early", "age@constant"))
})


# .hzr_stepwise_forward_step ------------------------------------------------

test_that("forward step enters the strongest candidate under Wald", {
  obj <- .fit_stepwise_base()
  step <- .hzr_stepwise_forward_step(
    obj$fit, scope = ~ x1 + x2 + x3, data = obj$data,
    criterion = "wald", slentry = 0.30
  )
  expect_true(step$accepted)
  expect_identical(step$variable, "x1")
  expect_lt(step$score, 0.30)
  expect_lt(step$p_value, 0.05)  # signal should be highly significant
  # all_scores should list every candidate that was tried
  expect_true(nrow(step$all_scores) >= 3L)
})

test_that("forward step returns unchanged when no candidate clears threshold", {
  obj <- .fit_stepwise_base(signal_beta = 0)  # pure noise
  step <- .hzr_stepwise_forward_step(
    obj$fit, scope = ~ x1 + x2 + x3, data = obj$data,
    criterion = "wald", slentry = 0.01   # strict threshold
  )
  expect_false(step$accepted)
  expect_identical(step$fit, obj$fit)
})

test_that("forward step under AIC uses ΔAIC rule and accepts when negative", {
  obj <- .fit_stepwise_base()
  step <- .hzr_stepwise_forward_step(
    obj$fit, scope = ~ x1 + x2 + x3, data = obj$data,
    criterion = "aic"
  )
  expect_true(step$accepted)
  expect_identical(step$variable, "x1")
  expect_lt(step$score, 0)        # ΔAIC improvement
  expect_equal(step$score, step$delta_aic)
})

test_that("forward step respects force_out", {
  obj <- .fit_stepwise_base()
  step <- .hzr_stepwise_forward_step(
    obj$fit, scope = NULL, data = obj$data,
    criterion = "wald", slentry = 0.30,
    force_out = "x1"   # the signal is locked out
  )
  expect_true(step$accepted || !step$accepted)  # may or may not accept,
  # but x1 must not be the chosen variable:
  expect_false(isTRUE(step$variable == "x1"))
})

test_that("forward step reports no candidates when scope is empty", {
  obj <- .fit_stepwise_base()
  step <- .hzr_stepwise_forward_step(
    obj$fit, scope = character(), data = obj$data,
    criterion = "wald"
  )
  expect_false(step$accepted)
  expect_identical(nrow(step$all_scores), 0L)
})
