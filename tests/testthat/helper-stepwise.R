# Shared fixtures for Wald / candidate-score / stepwise tests.  testthat
# auto-sources any helper-*.R file before running test-*.R files.

.fit_weibull <- function(n = 200L, betas = 0.3, seed = 42L) {
  set.seed(seed)
  X <- matrix(rnorm(n * length(betas)), ncol = length(betas))
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  lp <- as.numeric(X %*% betas)
  time <- rexp(n, rate = exp(lp))
  status <- rep(1L, n)
  hazard(
    time   = time,
    status = status,
    x      = X,
    theta  = c(0.5, 1.0, rep(0, ncol(X))),
    dist   = "weibull",
    fit    = TRUE
  )
}

# Used by the driver and integration tests: intercept-only Weibull
# base fit with a strong signal (x1) plus three noise covariates in
# the data frame, ready for forward selection.
.fit_driver_base <- function(n = 500L, seed = 911L, signal_beta = 0.9) {
  set.seed(seed)
  df <- data.frame(
    time   = rexp(n, rate = 1),
    status = rep(1L, n),
    x1     = rnorm(n),
    x2     = rnorm(n),
    x3     = rnorm(n),
    x4     = rnorm(n)
  )
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

# Overfitted Weibull including three covariates (one signal, two
# noise) — the starting point for backward-selection tests.
.fit_overfitted <- function(n = 400L, seed = 101L, signal_beta = 0.8) {
  set.seed(seed)
  df <- data.frame(
    time   = rexp(n, rate = 1),
    status = rep(1L, n),
    x1     = rnorm(n),
    x2     = rnorm(n),
    x3     = rnorm(n)
  )
  df$time <- df$time * exp(-signal_beta * df$x1)

  fit <- hazard(
    Surv(time, status) ~ x1 + x2 + x3,
    data = df,
    theta = c(0.5, 1.0, 0, 0, 0),
    dist = "weibull",
    fit = TRUE
  )
  list(fit = fit, data = df)
}

# Intercept-only Weibull used by the forward-step unit tests.
.fit_stepwise_base <- function(n = 400L, seed = 101L, signal_beta = 0.8) {
  set.seed(seed)
  df <- data.frame(
    time   = rexp(n, rate = 1),
    status = rep(1L, n),
    x1     = rnorm(n),
    x2     = rnorm(n),
    x3     = rnorm(n)
  )
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
