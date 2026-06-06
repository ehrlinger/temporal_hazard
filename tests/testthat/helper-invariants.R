# Tier-1 correctness invariants (see inst/dev/CORRECTNESS-STRATEGY.md).
#
# These are R-only, single-engine checks: properties every fitted model must
# satisfy by the mathematics, with no SAS reference and no captured fixtures.
# They are the cheap workhorse of the correctness regime -- they caught (or
# would have caught) the multiphase vcov, CoE left-truncation, deciles-method,
# and CoE final-adjustment bugs.
#
# Covered here: distributional sanity (S in [0,1] and monotone; H monotone and
# >= 0; S = exp(-H)), predict-type internal consistency, multiphase log-
# likelihood self-consistency (LL at returned theta == reported objective), and
# the Conservation-of-Events identity. Future expansion (documented, not yet
# implemented): stationarity (score ~ 0 at the MLE) and analytic-vs-numDeriv
# gradient agreement -- the latter is already exercised per-distribution by the
# analytic-Hessian tests.

# ---------------------------------------------------------------------------
# Synthetic data generators spanning the feature space.
# ---------------------------------------------------------------------------

# Right-censored survival data with one covariate. Deterministic given `seed`.
.inv_sim <- function(n = 300, rate = 0.3, cens_max = 6, seed = 1) {
  set.seed(seed)
  x <- rnorm(n)
  t <- rexp(n, rate = rate) + 0.01
  cens <- runif(n, 0.2, cens_max)
  data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x)
}

# Left-truncated (counting-process) two-phase-ish data. Deterministic.
.inv_sim_lt <- function(n = 400, seed = 7) {
  set.seed(seed)
  u <- runif(n)
  tev <- ifelse(u < 0.5, rexp(n, 3), rexp(n, 0.4))
  cens <- runif(n, 0.2, 6)
  stop <- pmin(tev, cens)
  data.frame(start = pmin(runif(n, 0, 1), 0.8 * stop),
             stop = stop, event = as.integer(tev <= cens))
}

# ---------------------------------------------------------------------------
# Model matrix: a list of fitted models spanning distributions, covariates,
# conservation, left-truncation and weights. Each entry carries the newdata
# grid (a fixed covariate profile over a time grid) for the distributional
# checks, and the raw data for the CoE / self-consistency checks.
# ---------------------------------------------------------------------------
.inv_models <- function() {
  df  <- .inv_sim()
  dlt <- .inv_sim_lt()
  set.seed(11)
  df$w <- sample(1:3, nrow(df), replace = TRUE)  # integer weights
  tg  <- seq(0.05, 5, length.out = 40)

  mk <- function(name, fit, grid, data = NULL, multiphase = FALSE,
                 conserve = FALSE, time = NULL, status = NULL,
                 time_lower = NULL, weights = NULL) {
    list(name = name, fit = fit, grid = grid, data = data,
         multiphase = multiphase, conserve = conserve,
         time = time, status = status, time_lower = time_lower,
         weights = weights)
  }

  models <- list()

  # --- Single-distribution families (intercept + one covariate) ------------
  models[[length(models) + 1L]] <- mk("weibull",
    hazard(survival::Surv(time, status) ~ x, data = df,
           theta = c(0.3, 1, 0), dist = "weibull", fit = TRUE),
    grid = data.frame(time = tg, x = 0))
  models[[length(models) + 1L]] <- mk("exponential",
    hazard(survival::Surv(time, status) ~ x, data = df,
           theta = c(log(0.3), 0), dist = "exponential", fit = TRUE),
    grid = data.frame(time = tg, x = 0))
  models[[length(models) + 1L]] <- mk("loglogistic",
    hazard(survival::Surv(time, status) ~ x, data = df,
           theta = c(0, 0.2, 0), dist = "loglogistic", fit = TRUE),
    grid = data.frame(time = tg, x = 0))
  models[[length(models) + 1L]] <- mk("lognormal",
    hazard(survival::Surv(time, status) ~ x, data = df,
           theta = c(1, 0, 0), dist = "lognormal", fit = TRUE),
    grid = data.frame(time = tg, x = 0))

  # --- Multiphase: conserve on/off, covariates, truncation, weights --------
  two_phase <- function(formula, data, conserve) {
    set.seed(101)
    hazard(formula, data = data, dist = "multiphase",
           phases = list(
             early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                  fixed = "shapes"),
             constant = hzr_phase("constant")),
           fit = TRUE,
           control = list(n_starts = 1, maxit = 500, conserve = conserve))
  }

  models[[length(models) + 1L]] <- mk("multiphase-coe",
    two_phase(survival::Surv(time, status) ~ 1, df, TRUE),
    grid = data.frame(time = tg), data = df, multiphase = TRUE, conserve = TRUE,
    time = df$time, status = df$status)

  models[[length(models) + 1L]] <- mk("multiphase-free",
    two_phase(survival::Surv(time, status) ~ 1, df, FALSE),
    grid = data.frame(time = tg), data = df, multiphase = TRUE, conserve = FALSE,
    time = df$time, status = df$status)

  set.seed(101)
  mp_cov <- hazard(survival::Surv(time, status) ~ 1, data = df,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ x),
      constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1, maxit = 500, conserve = TRUE))
  models[[length(models) + 1L]] <- mk("multiphase-covariate",
    mp_cov, grid = data.frame(time = tg, x = 0), data = df,
    multiphase = TRUE, conserve = TRUE, time = df$time, status = df$status)

  models[[length(models) + 1L]] <- mk("multiphase-truncated",
    two_phase(survival::Surv(start, stop, event) ~ 1, dlt, TRUE),
    grid = data.frame(time = tg), data = dlt, multiphase = TRUE,
    conserve = TRUE, time = dlt$stop, status = dlt$event,
    time_lower = dlt$start)

  set.seed(101)
  mp_w <- hazard(survival::Surv(time, status) ~ 1, data = df,
    dist = "multiphase", weights = df$w,
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1, maxit = 500, conserve = TRUE))
  models[[length(models) + 1L]] <- mk("multiphase-weighted",
    mp_w, grid = data.frame(time = tg), data = df, multiphase = TRUE,
    conserve = TRUE, time = df$time, status = df$status, weights = df$w)

  models
}

# ---------------------------------------------------------------------------
# Invariant assertions.
# ---------------------------------------------------------------------------

# S(t) in [0,1] and non-increasing; H(t) >= 0 and non-decreasing; S = exp(-H)
# (which also pins predict("survival") and predict("cumulative_hazard") to be
# mutually consistent).
.inv_assert_distribution <- function(m) {
  S <- unname(predict(m$fit, newdata = m$grid, type = "survival"))
  H <- unname(predict(m$fit, newdata = m$grid, type = "cumulative_hazard"))
  testthat::expect_true(all(is.finite(S)) && all(is.finite(H)),
                        label = paste0(m$name, ": finite S/H"))
  testthat::expect_true(all(S >= -1e-9 & S <= 1 + 1e-9),
                        label = paste0(m$name, ": S in [0,1]"))
  testthat::expect_true(all(diff(S) <= 1e-9),
                        label = paste0(m$name, ": S non-increasing"))
  testthat::expect_true(all(H >= -1e-9),
                        label = paste0(m$name, ": H >= 0"))
  testthat::expect_true(all(diff(H) >= -1e-9),
                        label = paste0(m$name, ": H non-decreasing"))
  testthat::expect_equal(S, exp(-H), tolerance = 1e-6,
                         label = paste0(m$name, ": S = exp(-H)"))
}

# LL recomputed at the returned coefficients equals the reported objective.
.inv_assert_self_consistent <- function(m) {
  w <- if (is.null(m$weights)) rep(1, length(m$time)) else m$weights
  ll <- TemporalHazard:::.hzr_logl_multiphase(
    m$fit$fit$theta, m$time, m$status, time_lower = m$time_lower,
    phases = m$fit$fit$phases, covariate_counts = m$fit$fit$covariate_counts,
    x_list = m$fit$fit$x_list, weights = w)
  testthat::expect_equal(ll, m$fit$fit$objective, tolerance = 1e-6,
                         label = paste0(m$name, ": LL(theta) == objective"))
}

# Conservation of Events: sum of (weighted) observed events equals sum of
# (weighted) predicted cumulative hazard over follow-up, H(stop) - H(start).
.inv_assert_coe <- function(m) {
  th <- m$fit$fit$theta
  ph <- m$fit$fit$phases
  cc <- m$fit$fit$covariate_counts
  xl <- m$fit$fit$x_list
  w  <- if (is.null(m$weights)) rep(1, length(m$time)) else m$weights
  H_stop  <- TemporalHazard:::.hzr_multiphase_cumhaz(m$time, th, ph, cc, xl)
  H_start <- if (is.null(m$time_lower)) 0 else
    TemporalHazard:::.hzr_multiphase_cumhaz(m$time_lower, th, ph, cc, xl)
  pred <- sum(w * (H_stop - H_start))
  obs  <- sum(w * m$status)
  testthat::expect_equal(pred, obs, tolerance = 1e-3 * obs,
                         label = paste0(m$name, ": CoE sum H == sum events"))
}
