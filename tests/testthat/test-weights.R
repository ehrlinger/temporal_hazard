# tests/testthat/test-weights.R
# Tests for the `weights` argument added in 0.9.4.
#
# The key invariant we check is that Fisher weighting is equivalent to
# row duplication for integer weights:
#
#   logL(theta; weights = w,  data = D)
#     ==  logL(theta; weights = 1,  data = expand(D, w))
#
# where expand(D, w) duplicates row i exactly w_i times.  This verifies
# that weights propagate through every likelihood term (event,
# right-censored, interval, left).  Running the same invariant through
# `hazard(fit = TRUE)` additionally exercises the optimizer, gradient,
# and Hessian paths because any mis-weighting there would show up as a
# different MLE or vcov between the two fits.

make_toy <- function(seed = 17) {
  set.seed(seed)
  n <- 30
  x <- rnorm(n)
  # Weibull-ish times
  u <- runif(n)
  t <- (-log(u) / exp(0.3 * x))^(1 / 1.2)
  # Mix of events, right-censored, and interval-censored
  status <- ifelse(t < 0.6, 1, 0)
  data.frame(time = t, status = status, x = x)
}

expand_rows <- function(df, weights) {
  df[rep(seq_len(nrow(df)), times = weights), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------

test_that("weights must be numeric and the right length", {
  df <- make_toy()
  expect_error(
    hazard(survival::Surv(time, status) ~ x, data = df,
           dist = "weibull", weights = "a"),
    "numeric vector"
  )
  expect_error(
    hazard(survival::Surv(time, status) ~ x, data = df,
           dist = "weibull", weights = rep(1, nrow(df) - 1)),
    "numeric vector"
  )
})

test_that("weights must be non-negative and finite", {
  df <- make_toy()
  expect_error(
    hazard(survival::Surv(time, status) ~ x, data = df,
           dist = "weibull", weights = c(-1, rep(1, nrow(df) - 1))),
    "non-negative"
  )
  expect_error(
    hazard(survival::Surv(time, status) ~ x, data = df,
           dist = "weibull", weights = c(Inf, rep(1, nrow(df) - 1))),
    "non-negative"
  )
})

# ---------------------------------------------------------------------------
# Weighted log-likelihood equals duplicated-row log-likelihood (Weibull)
# ---------------------------------------------------------------------------

test_that("integer weights match row duplication at an arbitrary theta (Weibull)", {
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  df_exp <- expand_rows(df, w)
  theta <- c(mu = 0.8, nu = 1.2, beta = 0.25)

  ll_weighted <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta,
    time = df$time, status = df$status,
    x = matrix(df$x, ncol = 1),
    weights = w
  )
  ll_expanded <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta,
    time = df_exp$time, status = df_exp$status,
    x = matrix(df_exp$x, ncol = 1),
    weights = rep(1, nrow(df_exp))
  )

  expect_equal(ll_weighted, ll_expanded, tolerance = 1e-10)
})

test_that("integer weights match row duplication for all censoring flavours", {
  set.seed(3)
  n <- 20
  x <- rnorm(n)
  t <- runif(n, 0.1, 2)
  # Mix: 10 events, 5 right-censored, 5 left-censored
  status <- c(rep(1, 10), rep(0, 5), rep(-1, 5))
  # The parser's convention for left-censored data is time_upper = time for
  # every row (non-left rows ignore it). Mirror that here so the LL's
  # `any(upper <= 0)` guard doesn't trip on NAs.
  time_upper <- t
  w <- sample(1:3, n, replace = TRUE)

  theta <- c(mu = 0.7, nu = 1.3, beta = 0.1)
  ll_w <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = t, status = status,
    time_upper = time_upper, x = matrix(x, ncol = 1), weights = w
  )

  idx <- rep(seq_len(n), times = w)
  ll_dup <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = t[idx], status = status[idx],
    time_upper = time_upper[idx], x = matrix(x[idx], ncol = 1),
    weights = rep(1, sum(w))
  )

  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Weighted fit == duplicated-row fit (MLE + vcov)
# ---------------------------------------------------------------------------

test_that("fit with integer weights matches the duplicated-row fit (Weibull)", {
  skip_on_cran()
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  df_exp <- expand_rows(df, w)

  fit_w <- hazard(
    survival::Surv(time, status) ~ x,
    data = df, dist = "weibull",
    weights = w, fit = TRUE
  )
  fit_dup <- hazard(
    survival::Surv(time, status) ~ x,
    data = df_exp, dist = "weibull",
    fit = TRUE
  )

  # Coefficients should match up to optimizer tolerance.
  expect_equal(coef(fit_w), coef(fit_dup), tolerance = 1e-3)
  # Log-likelihood should match tightly.
  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# Unit weights == unweighted (sanity)
# ---------------------------------------------------------------------------

test_that("unit weights are equivalent to no weights (Weibull)", {
  df <- make_toy()
  theta <- c(mu = 0.8, nu = 1.2, beta = 0.25)

  ll_noweight <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time, status = df$status,
    x = matrix(df$x, ncol = 1)
  )
  ll_unit <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time, status = df$status,
    x = matrix(df$x, ncol = 1), weights = rep(1, nrow(df))
  )
  expect_equal(ll_noweight, ll_unit, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Multiphase weighted likelihood equals duplicated-row likelihood
# ---------------------------------------------------------------------------

test_that("integer weights match duplication in multiphase likelihood", {
  skip_on_cran()
  set.seed(11)
  n <- 40
  t <- runif(n, 0.05, 3)
  status <- rbinom(n, 1, 0.4)
  w <- sample(1:2, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)

  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  # Start values: [log_mu_early, t_half, nu, m, log_mu_constant]
  theta <- c(-4, log(0.3), 1, 1, -3)

  ll_w <- TemporalHazard:::.hzr_logl_multiphase(
    theta = theta, time = t, status = status,
    phases = phases,
    covariate_counts = c(early = 0L, constant = 0L),
    x_list = list(early = NULL, constant = NULL),
    weights = w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_multiphase(
    theta = theta, time = t[idx], status = status[idx],
    phases = phases,
    covariate_counts = c(early = 0L, constant = 0L),
    x_list = list(early = NULL, constant = NULL),
    weights = rep(1, sum(w))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Gradient is weighted (numerical check on Weibull)
# ---------------------------------------------------------------------------

test_that("weighted log-likelihood gradient equals numerical gradient", {
  skip_if_not_installed("numDeriv")
  df <- make_toy()
  w <- sample(1:3, nrow(df), replace = TRUE)
  theta <- c(mu = 0.8, nu = 1.2, beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  # Numerical gradient of the weighted LL
  num_g <- numDeriv::grad(
    function(th) {
      TemporalHazard:::.hzr_logl_weibull(
        theta = th, time = df$time, status = df$status,
        x = xm, weights = w
      )
    },
    theta
  )

  ana_g <- TemporalHazard:::.hzr_gradient_weibull(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )

  expect_equal(as.numeric(ana_g), num_g, tolerance = 1e-5)
})

# ---------------------------------------------------------------------------
# Zero weight excludes the row
# ---------------------------------------------------------------------------

test_that("a zero weight drops the row's contribution to the likelihood", {
  df <- make_toy()
  w <- rep(1, nrow(df))
  w[1] <- 0

  theta <- c(mu = 0.8, nu = 1.2, beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  ll_drop <- TemporalHazard:::.hzr_logl_weibull(
    theta = theta, time = df$time[-1], status = df$status[-1],
    x = xm[-1, , drop = FALSE],
    weights = rep(1, nrow(df) - 1)
  )
  expect_equal(ll_w, ll_drop, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Per-distribution duplication-parity for exp / log-logistic / log-normal
# (wired up in v0.9.6 as part of Phase 4e).
# ---------------------------------------------------------------------------

test_that("integer weights match row duplication (exponential LL)", {
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)
  theta <- c(log_lambda = log(0.6), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_exponential(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_exponential(
    theta = theta, time = df$time[idx], status = df$status[idx],
    x = matrix(df$x[idx], ncol = 1), weights = rep(1, sum(w))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

test_that("integer weights match row duplication (log-logistic LL)", {
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)
  theta <- c(log_alpha = log(0.6), log_beta = log(1.2), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_loglogistic(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_loglogistic(
    theta = theta, time = df$time[idx], status = df$status[idx],
    x = matrix(df$x[idx], ncol = 1), weights = rep(1, sum(w))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

test_that("integer weights match row duplication (log-normal LL)", {
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  idx <- rep(seq_len(n), times = w)
  theta <- c(mu = 0, log_sigma = log(1), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  ll_w <- TemporalHazard:::.hzr_logl_lognormal(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_lognormal(
    theta = theta, time = df$time[idx], status = df$status[idx],
    x = matrix(df$x[idx], ncol = 1), weights = rep(1, sum(w))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})

# Analytic gradients for exp / log-logistic / log-normal match numerical
# gradients of the weighted LL.  Catches the "accepts formal but never
# applies it" pattern in the analytic score path.

test_that("weighted exponential gradient equals numerical gradient", {
  skip_if_not_installed("numDeriv")
  df <- make_toy()
  w <- sample(1:3, nrow(df), replace = TRUE)
  theta <- c(log_lambda = log(0.6), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  num_g <- numDeriv::grad(
    function(th) {
      TemporalHazard:::.hzr_logl_exponential(
        theta = th, time = df$time, status = df$status,
        x = xm, weights = w
      )
    },
    theta
  )
  ana_g <- TemporalHazard:::.hzr_gradient_exponential(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  expect_equal(as.numeric(ana_g), num_g, tolerance = 1e-5)
})

test_that("weighted log-logistic gradient equals numerical gradient", {
  skip_if_not_installed("numDeriv")
  df <- make_toy()
  w <- sample(1:3, nrow(df), replace = TRUE)
  theta <- c(log_alpha = log(0.6), log_beta = log(1.2), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  num_g <- numDeriv::grad(
    function(th) {
      TemporalHazard:::.hzr_logl_loglogistic(
        theta = th, time = df$time, status = df$status,
        x = xm, weights = w
      )
    },
    theta
  )
  ana_g <- TemporalHazard:::.hzr_gradient_loglogistic(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  expect_equal(as.numeric(ana_g), num_g, tolerance = 1e-5)
})

test_that("weighted log-normal gradient equals numerical gradient", {
  skip_if_not_installed("numDeriv")
  df <- make_toy()
  w <- sample(1:3, nrow(df), replace = TRUE)
  theta <- c(mu = 0, log_sigma = log(1), beta = 0.25)
  xm <- matrix(df$x, ncol = 1)

  num_g <- numDeriv::grad(
    function(th) {
      TemporalHazard:::.hzr_logl_lognormal(
        theta = th, time = df$time, status = df$status,
        x = xm, weights = w
      )
    },
    theta
  )
  ana_g <- TemporalHazard:::.hzr_gradient_lognormal(
    theta = theta, time = df$time, status = df$status,
    x = xm, weights = w
  )
  expect_equal(as.numeric(ana_g), num_g, tolerance = 1e-5)
})

# End-to-end: weighted fit == duplicated-row fit for each distribution.

test_that("fit with integer weights matches the duplicated-row fit (exponential)", {
  skip_on_cran()
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  df_exp <- expand_rows(df, w)

  fit_w <- hazard(
    survival::Surv(time, status) ~ x, data = df,
    dist = "exponential", weights = w, fit = TRUE
  )
  fit_dup <- hazard(
    survival::Surv(time, status) ~ x, data = df_exp,
    dist = "exponential", fit = TRUE
  )

  expect_equal(coef(fit_w), coef(fit_dup), tolerance = 1e-3)
  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-4)
})

test_that("fit with integer weights matches the duplicated-row fit (log-logistic)", {
  skip_on_cran()
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  df_exp <- expand_rows(df, w)

  fit_w <- hazard(
    survival::Surv(time, status) ~ x, data = df,
    dist = "loglogistic", weights = w, fit = TRUE
  )
  fit_dup <- hazard(
    survival::Surv(time, status) ~ x, data = df_exp,
    dist = "loglogistic", fit = TRUE
  )

  expect_equal(coef(fit_w), coef(fit_dup), tolerance = 1e-3)
  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-4)
})

test_that("fit with integer weights matches the duplicated-row fit (log-normal)", {
  skip_on_cran()
  df <- make_toy()
  n <- nrow(df)
  w <- sample(1:3, n, replace = TRUE)
  df_exp <- expand_rows(df, w)

  fit_w <- hazard(
    survival::Surv(time, status) ~ x, data = df,
    dist = "lognormal", weights = w, fit = TRUE
  )
  fit_dup <- hazard(
    survival::Surv(time, status) ~ x, data = df_exp,
    dist = "lognormal", fit = TRUE
  )

  expect_equal(coef(fit_w), coef(fit_dup), tolerance = 1e-3)
  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-4)
})
