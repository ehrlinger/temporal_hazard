# tests/testthat/test-phase-specific-covariates.R
#
# Phase 7d: correctness tests for phase-specific covariates (hzr_phase(formula=)).
#
# hzr_phase(formula = ~ ...) lets a covariate modulate the scale of a single
# phase, so the same covariate can carry different effects in different phases.
# The shipped smoke test only checks that such a fit runs; it does not check
# that the coefficients are correct or that a phase-specific covariate stays
# confined to its phase.
#
# These are simulation-based recovery tests -- the honest substitute for a SAS
# parity fixture when no reference output exists. They guard against this
# package's recurring "accepts the formal but never applies it" failure mode
# (seen with weights, counting-process times, and CoE weights) by confirming a
# covariate entered into a specific phase (a) recovers its true coefficient and
# (b) affects only that phase.

# ---------------------------------------------------------------------------
# Simulation helper: additive two-phase model (CDF early + constant), with
# optional phase-specific covariate effects on each phase's scale.
#   H(t | x) = mu_e(x) * Phi_early(t) + mu_c(x) * t
# Event times are drawn by inverting the additive cumulative hazard on a grid.
# ---------------------------------------------------------------------------
.sim_phase_cov <- function(n, beta_early = NULL, beta_const = NULL,
                           log_mu_e = -1.5, log_mu_c = -2.5,
                           t_half = 0.5, tmax = 80, seed = 1) {
  set.seed(seed)
  x      <- rnorm(n)
  tgrid  <- seq(0.001, tmax, length.out = 4000)
  Phi_e  <- hzr_phase_cumhaz(tgrid, t_half = t_half, nu = 1, m = 1, type = "cdf")
  U      <- -log(runif(n))
  sim_t  <- numeric(n)
  status <- integer(n)
  for (i in seq_len(n)) {
    mu_e <- exp(log_mu_e + if (!is.null(beta_early)) beta_early * x[i] else 0)
    mu_c <- exp(log_mu_c + if (!is.null(beta_const)) beta_const * x[i] else 0)
    H    <- mu_e * Phi_e + mu_c * tgrid
    # First grid index where the accumulated hazard reaches the target U[i];
    # NA only when the event never occurs within the grid (right-censored).
    j    <- match(TRUE, H >= U[i])
    if (is.na(j)) {
      sim_t[i]  <- max(tgrid)
      status[i] <- 0L            # right-censored: event beyond the grid
    } else {
      sim_t[i]  <- tgrid[j]
      status[i] <- 1L            # event observed (even if on the last point)
    }
  }
  data.frame(
    int_dead = sim_t,
    dead     = status,
    x        = x
  )
}

# ---------------------------------------------------------------------------
# Structural tests (fast, always run): a phase-specific covariate produces a
# coefficient for its own phase only.
# ---------------------------------------------------------------------------

test_that("a covariate in one phase produces a beta for that phase only", {
  df <- .sim_phase_cov(n = 200, beta_early = 0.8, seed = 3)

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = df,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant")
    ),
    fit     = TRUE,
    control = list(n_starts = 1L, maxit = 50L)
  )

  nm <- names(coef(fit))
  expect_true("early.x" %in% nm)
  expect_false("constant.x" %in% nm)
  # covariate_counts reflects the per-phase covariate placement
  expect_equal(unname(fit$fit$covariate_counts[["early"]]), 1L)
  expect_equal(unname(fit$fit$covariate_counts[["constant"]]), 0L)
})

test_that("the same covariate can enter two phases independently", {
  df <- .sim_phase_cov(n = 200, beta_early = 0.8, beta_const = -0.4, seed = 4)

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = df,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant", formula = ~ x)
    ),
    fit     = TRUE,
    control = list(n_starts = 1L, maxit = 50L)
  )

  nm <- names(coef(fit))
  expect_true(all(c("early.x", "constant.x") %in% nm))
  expect_equal(unname(fit$fit$covariate_counts[["early"]]), 1L)
  expect_equal(unname(fit$fit$covariate_counts[["constant"]]), 1L)
})

# ---------------------------------------------------------------------------
# Recovery tests (heavier simulations; skip on CRAN): the fitted coefficients
# match the data-generating truth.
# ---------------------------------------------------------------------------

test_that("phase-specific covariate in the early phase is recovered", {
  skip_on_cran()
  beta_early <- 0.8
  log_mu_e   <- -1.5
  log_mu_c   <- -3.0
  df <- .sim_phase_cov(n = 3000, beta_early = beta_early,
                       log_mu_e = log_mu_e, log_mu_c = log_mu_c, seed = 7)

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = df,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant")
    ),
    fit     = TRUE,
    control = list(n_starts = 3L, maxit = 1000L)
  )
  co <- coef(fit)

  expect_equal(unname(co[["early.x"]]),        beta_early, tolerance = 0.12)
  expect_equal(unname(co[["early.log_mu"]]),   log_mu_e,   tolerance = 0.15)
  expect_equal(unname(co[["constant.log_mu"]]), log_mu_c,  tolerance = 0.20)
})

test_that("the same covariate recovers opposite-sign effects across phases", {
  skip_on_cran()
  beta_early <- 1.0
  beta_const <- -0.6
  log_mu_e   <- -1.5
  log_mu_c   <- -2.5
  df <- .sim_phase_cov(n = 5000, beta_early = beta_early, beta_const = beta_const,
                       log_mu_e = log_mu_e, log_mu_c = log_mu_c, seed = 11)

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = df,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant", formula = ~ x)
    ),
    fit     = TRUE,
    control = list(n_starts = 3L, maxit = 1000L)
  )
  co <- coef(fit)

  # The defining feature: one covariate, two different (here opposite-sign)
  # phase-specific effects, both recovered.
  expect_equal(unname(co[["early.x"]]),    beta_early, tolerance = 0.12)
  expect_equal(unname(co[["constant.x"]]), beta_const, tolerance = 0.12)
  expect_true(co[["early.x"]] > 0)
  expect_true(co[["constant.x"]] < 0)
})

test_that("a covariate confined to one phase does not leak into the other", {
  skip_on_cran()
  # Truth: x affects ONLY the early phase. A model that also offers x to the
  # constant phase should estimate the constant-phase effect near zero.
  df <- .sim_phase_cov(n = 5000, beta_early = 1.0, beta_const = NULL,
                       log_mu_e = -1.5, log_mu_c = -2.5, seed = 13)

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = df,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant", formula = ~ x)
    ),
    fit     = TRUE,
    control = list(n_starts = 3L, maxit = 1000L)
  )
  co <- coef(fit)

  expect_equal(unname(co[["early.x"]]),    1.0, tolerance = 0.15)
  expect_equal(unname(co[["constant.x"]]), 0.0, tolerance = 0.15)
})
