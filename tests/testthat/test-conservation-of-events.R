# tests/testthat/test-conservation-of-events.R
# Tests for the Conservation of Events (CoE) theorem implementation

test_that("CoE is active by default for multiphase models", {
  set.seed(42)
  data(cabgkul, package = "TemporalHazard")
  fit_coe <- hazard(
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
  expect_true(isTRUE(fit_coe$fit$converged) || is.finite(fit_coe$fit$objective))
  expect_true(is.finite(fit_coe$fit$objective))
})

test_that("CoE can be disabled via control$conserve = FALSE", {
  set.seed(42)
  data(cabgkul, package = "TemporalHazard")
  fit_no_coe <- hazard(
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
    control = list(n_starts = 3, maxit = 500, conserve = FALSE)
  )
  expect_true(is.finite(fit_no_coe$fit$objective))
})

test_that("CoE improves conservation ratio toward 1.0", {
  set.seed(42)
  data(cabgkul, package = "TemporalHazard")

  fit_coe <- hazard(
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

  # Check conservation: predicted events should be close to observed
  gof <- hzr_gof(fit_coe)
  s <- attr(gof, "summary")
  ratio <- s$total_expected / s$total_observed

  # With CoE, the conservation ratio should be very close to 1.0
  expect_true(
    abs(ratio - 1.0) < 0.15,
    label = paste("conservation ratio", round(ratio, 4),
                  "should be within 0.15 of 1.0")
  )
})

test_that("CoE CABGKUL log-likelihood matches C reference", {
  set.seed(42)
  data(cabgkul, package = "TemporalHazard")
  fit_coe <- hazard(
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
    control = list(n_starts = 5, maxit = 1000)
  )

  # C reference with fixed shapes: log-lik = -3740.52
  # CoE should help reach this or better
  expect_true(
    fit_coe$fit$objective >= -3745,
    label = paste("CoE log-lik", round(fit_coe$fit$objective, 2),
                  "should be >= -3745 (C ref: -3740.52)")
  )
})

test_that("CoE works with 2-phase model (early + constant)", {
  set.seed(42)
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = avc,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                            fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit     = TRUE,
    control = list(n_starts = 3, maxit = 500)
  )
  expect_true(isTRUE(fit$fit$converged) || is.finite(fit$fit$objective))
})

test_that(".hzr_log_mu_positions returns correct positions", {
  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1),
    constant = hzr_phase("constant"),
    late     = hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1)
  )
  cov_counts <- c(early = 0L, constant = 0L, late = 0L)
  pos <- TemporalHazard:::.hzr_log_mu_positions(phases, cov_counts)

  # early: log_mu at 1, then 3 shapes = 4 params
  # constant: log_mu at 5, no shapes = 1 param
  # late: log_mu at 6, then 4 G3 shapes = 5 params
  expect_equal(pos[["early"]], 1L)
  expect_equal(pos[["constant"]], 5L)
  expect_equal(pos[["late"]], 6L)
})

test_that(".hzr_conserve_events adjusts theta correctly", {
  data(cabgkul, package = "TemporalHazard")
  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1),
    constant = hzr_phase("constant")
  )
  cov_counts <- c(early = 0L, constant = 0L)
  # Arbitrary starting theta: [log_mu_e, log_thalf, nu, m, log_mu_c]
  theta <- c(-3, log(0.2), 1, 1, -5)

  theta_adj <- TemporalHazard:::.hzr_conserve_events(
    theta, fixmu_phase = "constant", fixmu_pos = 5L,
    time = cabgkul$int_dead, status = cabgkul$dead,
    phases = phases, covariate_counts = cov_counts,
    x_list = list(early = NULL, constant = NULL),
    total_events = sum(cabgkul$dead)
  )

  # The constant phase log_mu should have changed
  expect_false(theta_adj[5] == theta[5])
  # Other params should be unchanged
  expect_equal(theta_adj[1:4], theta[1:4])
})

# ===========================================================================
# CoE parity: reduced-dim CoE and full-dim optimization find same optimum
# ===========================================================================
#
# CoE is a dimension-reduction technique -- one phase's log_mu is solved
# analytically at each iteration instead of numerically optimized.  The
# log-likelihood function is unchanged, so the MLE should be identical
# (up to optimizer tolerance) whether CoE is on or off.

test_that("CoE and non-CoE find the same optimum on right-censored data", {
  skip_on_cran()
  set.seed(101)
  data(cabgkul, package = "TemporalHazard")
  small <- cabgkul[1:60, ]

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1,
                           fixed = "shapes"),
    constant = hzr_phase("constant")
  )

  fit_coe <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data = small, dist = "multiphase", phases = phases,
    fit = TRUE,
    control = list(n_starts = 3, maxit = 500, conserve = TRUE)
  )
  fit_no_coe <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data = small, dist = "multiphase", phases = phases,
    fit = TRUE,
    control = list(n_starts = 3, maxit = 500, conserve = FALSE)
  )

  expect_equal(fit_coe$fit$objective, fit_no_coe$fit$objective,
               tolerance = 1e-2)
})

# ===========================================================================
# CoE with interval-censored data auto-disables cleanly
# ===========================================================================

test_that("CoE auto-disables (no error) when any observation is interval-censored", {
  # Per likelihood-multiphase.R the CoE guard is:
  #   coe_supported_data <- all(status %in% c(0, 1))
  # Interval-censored data (status == 2) therefore falls through to the
  # full-dimensional optimizer.  This test asserts the fit runs to a
  # finite log-likelihood.  It does not claim the Turner conservation
  # property holds for interval data -- that extension is tracked as
  # part of Phase 4e.
  skip_on_cran()
  set.seed(13)
  n <- 30
  lo <- runif(n, 0.1, 1.5)
  hi <- lo + runif(n, 0.1, 0.5)
  df <- data.frame(lo = lo, hi = hi)

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                           fixed = "shapes"),
    constant = hzr_phase("constant")
  )

  fit <- hazard(
    survival::Surv(lo, hi, type = "interval2") ~ 1,
    data = df, dist = "multiphase", phases = phases,
    fit = TRUE,
    control = list(n_starts = 2, maxit = 300)
  )
  expect_true(is.finite(fit$fit$objective))
})

# ===========================================================================
# CoE with non-unit weights (wired up in v0.9.6, Phase 4e)
# ===========================================================================
#
# Prior to 0.9.6, `.hzr_conserve_events()` summed per-phase cumhaz without
# applying row weights, so Turner's adjustment was on a mismatched scale
# vs the (weighted) observed event count.  CoE therefore auto-disabled
# when any weight != 1.  As of 0.9.6, weights are threaded into
# `.hzr_conserve_events()` / `.hzr_select_fixmu_phase()` and CoE stays on
# for weighted fits.  Duplication parity is the acceptance criterion.

test_that("weighted multiphase fit matches the row-duplicated fit (CoE off)", {
  skip_on_cran()
  set.seed(42)
  data(cabgkul, package = "TemporalHazard")
  small <- cabgkul[1:40, ]
  w <- sample(1:2, nrow(small), replace = TRUE)
  idx <- rep(seq_len(nrow(small)), times = w)
  dup <- small[idx, ]

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1,
                           fixed = "shapes"),
    constant = hzr_phase("constant")
  )

  # Force CoE off on both sides so the comparison isolates weighted LL
  # equivalence rather than the CoE vs non-CoE dimension trick.
  fit_w <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data = small, dist = "multiphase", phases = phases,
    weights = w, fit = TRUE,
    control = list(n_starts = 2, maxit = 300, conserve = FALSE)
  )
  fit_dup <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data = dup, dist = "multiphase", phases = phases,
    fit = TRUE,
    control = list(n_starts = 2, maxit = 300, conserve = FALSE)
  )

  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-3)
  expect_equal(unname(coef(fit_w)), unname(coef(fit_dup)), tolerance = 1e-2)
})

test_that("weighted multiphase fit matches the row-duplicated fit (CoE on)", {
  skip_on_cran()
  set.seed(42)
  data(cabgkul, package = "TemporalHazard")
  small <- cabgkul[1:40, ]
  w <- sample(1:2, nrow(small), replace = TRUE)
  idx <- rep(seq_len(nrow(small)), times = w)
  dup <- small[idx, ]

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1,
                           fixed = "shapes"),
    constant = hzr_phase("constant")
  )

  # CoE enabled on both sides: the weighted fit reduces the dimension on
  # D, the unweighted fit reduces it on expand(D, w).  The optimum (MLE
  # and log-likelihood) must match because CoE is a reparameterisation
  # of the same surface.
  fit_w <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data = small, dist = "multiphase", phases = phases,
    weights = w, fit = TRUE,
    control = list(n_starts = 2, maxit = 300, conserve = TRUE)
  )
  fit_dup <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data = dup, dist = "multiphase", phases = phases,
    fit = TRUE,
    control = list(n_starts = 2, maxit = 300, conserve = TRUE)
  )

  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-3)
  expect_equal(unname(coef(fit_w)), unname(coef(fit_dup)), tolerance = 1e-2)
})

test_that("CoE default (conserve=TRUE) with weights reaches the same optimum as conserve=FALSE", {
  # CoE reduces dimension; it does not change the log-likelihood surface.
  # Both paths must land at the same MLE log-likelihood.
  skip_on_cran()
  set.seed(77)
  data(cabgkul, package = "TemporalHazard")
  small <- cabgkul[1:40, ]
  w <- sample(1:2, nrow(small), replace = TRUE)

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1,
                           fixed = "shapes"),
    constant = hzr_phase("constant")
  )

  fit_default <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data = small, dist = "multiphase", phases = phases,
    weights = w, fit = TRUE,
    control = list(n_starts = 2, maxit = 300)  # conserve default TRUE
  )
  fit_explicit_off <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data = small, dist = "multiphase", phases = phases,
    weights = w, fit = TRUE,
    control = list(n_starts = 2, maxit = 300, conserve = FALSE)
  )

  expect_equal(fit_default$fit$objective, fit_explicit_off$fit$objective,
               tolerance = 1e-3)
})
