# ---------------------------------------------------------------------------
# Conservation of Events under left truncation (counting-process entry times).
#
# The multiphase log-likelihood scores mu via Sum E = Sum [H(stop) - H(start)]
# for left-truncated rows (Surv(start, stop, event)). The Conservation-of-Events
# reparameterization must solve the SAME constraint, i.e. it must subtract the
# entry-time cumulative hazard H(start). If CoE instead conserves Sum H(stop)
# (ignoring the truncation term), it constrains the fit to the wrong manifold
# and the attained maximum log-likelihood drops below the unconstrained fit.
#
# Invariant exploited here: CoE only removes a redundant parameter direction, so
# objective(conserve = TRUE) must equal objective(conserve = FALSE). This holds
# trivially for untruncated data (start = 0) and is the discriminating check for
# truncated data. This is the synthetic analogue of the hz.te123.OMC fit-1 gap
# (gap-list P1 #6).
# ---------------------------------------------------------------------------

# Simulate a fittable two-phase data set with optional left truncation. The
# exact data-generating process is unimportant; what matters is that (a) both
# the CoE and free fits converge and (b) entry times are a non-negligible
# fraction of follow-up so the omitted H(start) term is material.
.sim_coe_lt <- function(n = 600, seed = 101, truncate = TRUE) {
  set.seed(seed)
  u <- runif(n)
  # Early burst of events plus a long constant-ish tail -> two phases.
  t_event <- ifelse(u < 0.5, rexp(n, rate = 3), rexp(n, rate = 0.4))
  cens    <- runif(n, 0.2, 6)
  stop    <- pmin(t_event, cens)
  event   <- as.integer(t_event <= cens)
  start   <- if (truncate) pmin(runif(n, 0, 1.0), 0.8 * stop) else rep(0, n)
  data.frame(start = start, stop = stop, event = event)
}

.coe_phases <- function() {
  list(
    early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
}

.coe_fit <- function(d, conserve, seed = 202) {
  # Reset the RNG before fitting so the multi-start perturbations are
  # deterministic: the CoE and free fits must not drift into different local
  # optima from back-to-back RNG draws (which would make the objective
  # comparison flaky across platforms/R versions).
  set.seed(seed)
  hazard(
    survival::Surv(start, stop, event) ~ 1,
    data    = d,
    dist    = "multiphase",
    phases  = .coe_phases(),
    fit     = TRUE,
    control = list(n_starts = 5, maxit = 1500, reltol = 1e-9, conserve = conserve)
  )
}

# Recompute the multiphase log-likelihood at a fitted object's returned
# coefficients (no further CoE applied). A self-consistent fit must reproduce
# its own reported objective here.
.coe_ll_at_par <- function(fit, d) {
  TemporalHazard:::.hzr_logl_multiphase(
    theta = fit$fit$theta, time = d$stop, status = d$event,
    time_lower = d$start, phases = fit$fit$phases,
    covariate_counts = fit$fit$covariate_counts, x_list = fit$fit$x_list
  )
}

test_that("CoE conserves events on the entry-time scale (left-truncated)", {
  skip_on_cran()
  skip_if_not_installed("survival")

  d <- .sim_coe_lt(truncate = TRUE)
  fit_coe  <- .coe_fit(d, conserve = TRUE)
  fit_free <- .coe_fit(d, conserve = FALSE)

  expect_true(isTRUE(fit_coe$fit$converged))
  expect_true(isTRUE(fit_free$fit$converged))

  # A correct CoE cannot lower the attainable maximum log-likelihood.
  expect_equal(fit_coe$fit$objective, fit_free$fit$objective,
               tolerance = 1e-3,
               label = "objective(conserve=TRUE) vs conserve=FALSE, left-truncated")

  # The returned coefficients must be self-consistent with the objective: the
  # final post-optimization CoE adjustment has to solve the conserved log_mu on
  # the same entry-time scale, else the returned mu drifts off the optimum.
  expect_equal(.coe_ll_at_par(fit_coe, d), fit_coe$fit$objective,
               tolerance = 1e-6,
               label = "LL at returned coefficients vs reported objective")
})

test_that("CoE matches the free fit for untruncated data (control)", {
  skip_on_cran()
  skip_if_not_installed("survival")

  d <- .sim_coe_lt(truncate = FALSE)
  fit_coe  <- .coe_fit(d, conserve = TRUE)
  fit_free <- .coe_fit(d, conserve = FALSE)

  expect_equal(fit_coe$fit$objective, fit_free$fit$objective,
               tolerance = 1e-3,
               label = "objective(conserve=TRUE) vs conserve=FALSE, untruncated")
})
