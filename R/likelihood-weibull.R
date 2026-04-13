#' @importFrom stats optim setNames
#' @keywords internal
NULL

# likelihood-weibull.R — Weibull parametric hazard likelihood, gradient, and optimizer
#
# MODEL
# -----
# Proportional-hazards (PH) parameterization:
#
#   h(t | x)  = μ ν t^(ν−1) exp(η)         hazard
#   H(t | x)  = (μ t)^ν exp(η)              cumulative hazard
#   S(t | x)  = exp(−H(t | x))              survival
#
# where μ > 0 (scale), ν > 0 (shape), η = x β (linear predictor).
#
# THETA LAYOUT
# ------------
#   theta[1] = μ   (must be > 0; enforced by L-BFGS-B lower bound)
#   theta[2] = ν   (must be > 0; enforced by L-BFGS-B lower bound)
#   theta[3:p] = β (covariate coefficients; unrestricted)
#
# Note: unlike the other distributions, Weibull optimises on the natural scale
# and uses L-BFGS-B with lower = 1e-6 instead of log-reparameterisation.
# This preserves direct interpretability of the μ/ν estimates.
#
# FUNCTIONS
# ---------
#   .hzr_logl_weibull()     — log-likelihood (optionally returning gradient)
#   .hzr_gradient_weibull() — analytical score vector
#   .hzr_optim_weibull()    — L-BFGS-B wrapper with post-fit Hessian SE

#' Weibull Parametric Hazard Likelihood
#'
#' Evaluate the log-likelihood and its derivatives for Weibull hazard models.
#' This is the core mathematical engine for model fitting.
#'
#' @keywords internal

#' Log-likelihood for Weibull hazard with covariates
#'
#' Computes the negative log-likelihood for a sample under the Weibull
#' parametric hazard model with optional linear-predictor covariates.
#'
#' @param theta Vector of parameters:
#'   theta\[1:n_shape\]: Shape/scale parameters (distribution-specific)
#'   theta\[n_shape+1:length\]: Covariate coefficients (linear on log-hazard scale)
#'   where n_shape = number of shape parameters (e.g., 2 for Weibull scale/shape)
#'
#' @param time Numeric vector of follow-up times (n)
#' @param status Numeric vector of event indicators: 1 = event, 0 = censored (n)
#' @param time_lower Optional numeric lower bound vector for interval-censored rows.
#'   Defaults to time if NULL.
#' @param time_upper Optional numeric upper bound vector for left/interval-censored rows.
#'   Defaults to time if NULL.
#' @param x Design matrix of covariates (n × p_coef); NULL for no covariates
#' @param dist_name Character name of baseline distribution ("weibull", etc.)
#' @param return_gradient Logical; if TRUE, attach gradient vector as attribute
#' @param return_hessian Logical; if TRUE, attach Hessian matrix as attribute
#'
#' @return Scalar log-likelihood value. If return_gradient = TRUE, gradient vector
#' is attached as attribute \code{"gradient"}. If return_hessian = TRUE, Hessian
#' matrix is attached as \code{"hessian"}.
#'
#' @details
#' The Weibull hazard model is parameterized as:
#'
#' \deqn{h(t | x) = \mu \nu t^{\nu-1} \exp(\eta)}
#'
#' where:
#' - \eqn{\mu > 0} is scale (MU)
#' - \eqn{\nu > 0} is shape (NU)
#' - \eqn{\eta = x \beta} is log-relative-hazard (covariates)
#'
#' The log-likelihood for right-censored data is:
#'
#' \deqn{\ell(\theta) = \sum_{i: \delta_i = 1} \log h(t_i | x_i)
#'   - \sum_i S(t_i | x_i)}
#'
#' where S is the survival function (1 - CDF).
#'
#' Mixed censoring status coding:
#' - 1: exact event at time
#' - 0: right-censored at time
#' - -1: left-censored with upper bound time_upper \(or time\)
#' - 2: interval-censored in the interval \(time_lower, time_upper\)
#'
#' @noRd
.hzr_logl_weibull <- function(
    theta,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    dist_name = "weibull",
    return_gradient = FALSE,
    return_hessian = FALSE) {

  # Shape parameter count for Weibull
  n_shape <- 2  # mu, nu (scale, shape)
  n <- length(time)

  # Extract parameters
  mu <- theta[1]    # scale > 0
  nu <- theta[2]    # shape > 0

  # Covariate coefficients (if any)
  if (!is.null(x)) {
    if (is.null(attr(x, "dimnames")[[2]])) {
      colnames(x) <- paste0("beta_", seq_len(ncol(x)))
    }
    beta <- theta[3:length(theta)]
    eta <- as.numeric(x %*% beta)  # n-vector
  } else {
    eta <- rep(0, n)
  }

  # Stability checks
  if (mu <= 0 || nu <= 0) {
    return(Inf)  # Infeasible: return large penalty
  }

  # Normalize censoring bounds to lower/upper vectors for unified formulas.
  lower <- if (is.null(time_lower)) time else time_lower
  upper <- if (is.null(time_upper)) time else time_upper

  # Exact events require strictly positive time because log(h(t)) is used.
  # Right-censored observations can be at time 0 (valid: log S(0) = 0).
  if (any(status == 1 & time <= 0)) {
    return(Inf)
  }
  if (any(status %in% c(-1, 2)) && any(upper <= 0)) {
    return(Inf)
  }
  if (any(status == 2) && any(lower[status == 2] >= upper[status == 2])) {
    return(Inf)
  }

  # Event-time hazard/cumulative-hazard (exact events only).
  haz_event <- mu * nu * (time ^ (nu - 1)) * exp(eta)
  cumhaz_event <- (mu * time) ^ nu * exp(eta)

  # Lower/upper cumulative hazards for censoring contributions.
  cumhaz_lower <- (mu * lower) ^ nu * exp(eta)
  cumhaz_upper <- (mu * upper) ^ nu * exp(eta)

  idx_event <- status == 1
  idx_right <- status == 0
  idx_left <- status == -1
  idx_interval <- status == 2

  # Exact events: log h(t) + log S(t) = log h(t) - H(t)
  ll_event <- if (any(idx_event)) {
    sum(log(haz_event[idx_event]) - cumhaz_event[idx_event])
  } else {
    0
  }

  # Right-censored: log S(t) = -H(t)
  ll_right <- if (any(idx_right)) {
    -sum(cumhaz_event[idx_right])
  } else {
    0
  }

  # Left-censored: log F(u) = log(1 - S(u)) = log(1 - exp(-H(u)))
  ll_left <- if (any(idx_left)) {
    sum(hzr_log1mexp(cumhaz_upper[idx_left]))
  } else {
    0
  }

  # Interval-censored: log(S(l) - S(u)) = -H(l) + log(1 - exp(-(H(u)-H(l))))
  ll_interval <- if (any(idx_interval)) {
    delta_h <- cumhaz_upper[idx_interval] - cumhaz_lower[idx_interval]
    sum(-cumhaz_lower[idx_interval] + hzr_log1mexp(delta_h))
  } else {
    0
  }

  logl <- ll_event + ll_right + ll_left + ll_interval

  if (!is.finite(logl)) {
    return(Inf)
  }

  # If gradient requested, compute score vector.
  if (return_gradient) {
    grad <- .hzr_gradient_weibull(
      theta, time, status, time_lower, time_upper, x, eta,
      cumhaz_event, haz_event, mu, nu, n_shape
    )
    attr(logl, "gradient") <- grad
  }

  # If Hessian requested (stub for M2)
  if (return_hessian) {
    attr(logl, "hessian") <- matrix(NA, length(theta), length(theta))
  }

  logl
}

#' Score vector (gradient of log-likelihood)
#'
#' Computes the score vector (gradient) of the Weibull log-likelihood w.r.t. all parameters.
#'
#' @param time_lower Optional lower bounds for interval-censored rows.
#' @param time_upper Optional upper bounds for left/interval-censored rows.
#'
#' The log-likelihood is:
#'   L = sum(delta_i * log h(t_i)) - sum(H(t_i))
#'
#' where h is hazard and H is cumulative hazard. Derivatives are:
#'
#' dL/dmu  = sum(delta_i / mu) - (nu / mu) * sum(H(t_i))
#' dL/dnu  = sum(delta_i / nu) + sum(delta_i * log(t_i)) - sum(log(mu * t_i) * H(t_i))
#' dL/dbeta_j = sum(delta_i * x_ij) - sum(H(t_i) * x_ij)  = t(X) %*% (delta - H)
#'
#' @noRd
.hzr_gradient_weibull <- function(
    theta,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    eta = NULL,
    cumhaz = NULL,
    haz = NULL,
    mu = NULL,
    nu = NULL,
    n_shape = 2) {

  n <- length(time)
  p <- length(theta)
  grad <- numeric(p)

  # Interval/left censoring currently uses a robust numerical gradient fallback.
  # The closed-form derivatives for these terms are deferred to a later pass.
  if (any(status %in% c(-1, 2))) {
    return(.hzr_numeric_grad_weibull(theta, time, status, time_lower, time_upper, x))
  }

  # Sanity check: if eta, cumhaz, haz not provided, compute them
  if (is.null(eta) || is.null(cumhaz) || is.null(haz)) {
    # Need to recompute from theta
    mu <- theta[1]
    nu <- theta[2]

    if (!is.null(x)) {
      beta <- theta[3:length(theta)]
      eta <- as.numeric(x %*% beta)
    } else {
      eta <- rep(0, n)
    }

    haz <- mu * nu * (time ^ (nu - 1)) * exp(eta)
    cumhaz <- (mu * time) ^ nu * exp(eta)
  }

  # ===== Gradient w.r.t. mu (scale parameter) =====
  # dL/dmu = sum(delta / mu) - (nu / mu) * sum(H)
  grad[1] <- sum(status) / mu - (nu / mu) * sum(cumhaz)

  # ===== Gradient w.r.t. nu (shape parameter) =====
  # dL/dnu = sum(delta / nu) + sum(delta * log(t)) - sum(log(mu*t) * H)
  grad[2] <- sum(status) / nu + sum(status * log(time)) - sum(log(mu * time) * cumhaz)

  # ===== Gradient w.r.t. beta (covariate coefficients) =====
  if (n_shape < p && !is.null(x)) {
    # Compute residual: delta - cumulative hazard
    residual <- status - cumhaz

    # dL/dbeta = t(X) %*% (delta - H)
    grad[3:p] <- as.numeric(crossprod(x, residual))
  }

  grad
}

#' @noRd
.hzr_optim_weibull <- function(
    time, status, time_lower = NULL, time_upper = NULL,
    x = NULL, theta_start, control = list()) {

  # ── Internal reparameterisation ────────────────────────────────────────────
  # The original hazard C code parameterises the cumulative hazard as
  #
  #     H(t | x) = exp(α + Xβ) · t^γ
  #
  # where α is a free intercept, γ = exp(ψ) is the shape, and β enters
  # the same log-hazard-rate space as α.  This decouples scale from shape
  # and avoids the degenerate ridge that appears when covariates with large
  # means (e.g. age ≈ 62) interact with the (log mu, log nu) parameterisation.
  #
  # We optimise over φ = (α, ψ, β₁, …) — all unconstrained — then
  # back-transform to the user-facing (μ, ν, β) scale:
  #
  #     ν = exp(ψ)      μ = exp(α / ν)
  #
  # and apply the delta method for the variance–covariance matrix.
  # ──────────────────────────────────────────────────────────────────────────

  n_shape <- 2L
  p <- length(theta_start)
  n <- length(time)

  # ── Convert user theta (mu, nu, beta) → internal (alpha, psi, beta) ──────
  mu_start    <- unname(theta_start[1])
  nu_start    <- unname(theta_start[2])
  alpha_start <- nu_start * log(mu_start)   # α = ν · log(μ)
  psi_start   <- log(nu_start)              # ψ = log(ν)
  beta_start  <- if (p > n_shape) unname(theta_start[(n_shape + 1):p]) else numeric(0)
  phi_start   <- c(alpha_start, psi_start, beta_start)

  # ── Internal likelihood: H(t|x) = exp(α + Xβ) · t^exp(ψ) ────────────────
  logl_internal <- function(theta, time, status, time_lower, time_upper, x, ...) {
    alpha <- theta[1]
    g     <- exp(theta[2])   # nu = shape
    beta  <- if (p > n_shape) theta[(n_shape + 1):p] else numeric(0)
    eta   <- if (!is.null(x) && length(beta) > 0) as.numeric(x %*% beta) else rep(0, n)

    # Left / interval censoring: delegate to full Weibull likelihood via
    # the natural-scale function, converting params on the fly.
    if (any(status %in% c(-1, 2))) {
      mu <- exp(alpha / g)
      return(.hzr_logl_weibull(c(mu, g, beta), time, status,
                               time_lower, time_upper, x, ...))
    }

    log_t <- log(pmax(time, .Machine$double.xmin))
    H     <- exp(alpha + eta) * (time ^ g)

    idx_event <- status == 1
    idx_right <- status == 0

    ll_event <- if (any(idx_event)) {
      sum(alpha + eta[idx_event] + theta[2] +
          (g - 1) * log_t[idx_event] - H[idx_event])
    } else 0

    ll_right <- if (any(idx_right)) -sum(H[idx_right]) else 0

    logl <- ll_event + ll_right
    if (!is.finite(logl)) return(Inf)
    logl
  }

  # ── Internal gradient ─────────────────────────────────────────────────────
  grad_internal <- function(theta, time, status, time_lower, time_upper, x, ...) {
    alpha <- theta[1]
    psi   <- theta[2]
    g     <- exp(psi)
    beta  <- if (p > n_shape) theta[(n_shape + 1):p] else numeric(0)
    eta   <- if (!is.null(x) && length(beta) > 0) as.numeric(x %*% beta) else rep(0, n)

    # Left / interval censoring → numerical gradient fallback.
    if (any(status %in% c(-1, 2))) {
      obj <- function(th) logl_internal(th, time, status, time_lower, time_upper, x)
      if (requireNamespace("numDeriv", quietly = TRUE))
        return(as.numeric(numDeriv::grad(obj, theta)))
      eps <- 1e-6
      grd <- numeric(p)
      for (i in seq_along(theta)) {
        tp <- tm <- theta
        tp[i] <- tp[i] + eps
        tm[i] <- tm[i] - eps
        grd[i] <- (obj(tp) - obj(tm)) / (2 * eps)
      }
      return(grd)
    }

    log_t <- log(pmax(time, .Machine$double.xmin))
    H     <- exp(alpha + eta) * (time ^ g)
    grad  <- numeric(p)

    # dL/dα = n_events − Σ H
    grad[1] <- sum(status) - sum(H)

    # dL/dψ = n_events + g · Σ_events log(t) − g · Σ H·log(t)
    grad[2] <- sum(status) + g * sum(status * log_t) - g * sum(H * log_t)

    # dL/dβ_j = X_j'(δ − H)
    if (p > n_shape && !is.null(x)) {
      grad[(n_shape + 1):p] <- as.numeric(crossprod(x, status - H))
    }

    grad
  }

  # ── Optimise (all unconstrained) ──────────────────────────────────────────
  result <- .hzr_optim_generic(
    logl_fn     = logl_internal,
    gradient_fn = grad_internal,
    time        = time, status = status,
    time_lower  = time_lower, time_upper = time_upper,
    x           = x,
    theta_start = phi_start,
    control     = control,
    use_bounds  = FALSE
  )

  # ── Back-transform to natural scale (μ, ν, β) ────────────────────────────
  alpha_hat <- result$par[1]
  psi_hat   <- result$par[2]
  nu_hat    <- exp(psi_hat)
  mu_hat    <- exp(alpha_hat / nu_hat)
  beta_hat  <- if (p > n_shape) result$par[(n_shape + 1):p] else numeric(0)

  result$par <- setNames(c(mu_hat, nu_hat, beta_hat), names(theta_start))

  # Delta method: Cov(μ, ν, β) = J · Cov(α, ψ, β) · Jᵀ
  #
  #   dμ/dα = μ / ν             dμ/dψ = −μ α / ν
  #   dν/dα = 0                 dν/dψ = ν
  #   dβ/dα = 0                 dβ/dψ = 0         dβ/dβ = I
  if (is.matrix(result$vcov) && !anyNA(result$vcov)) {
    J <- diag(p)
    J[1, 1] <- mu_hat / nu_hat
    J[1, 2] <- -mu_hat * alpha_hat / nu_hat
    J[2, 1] <- 0
    J[2, 2] <- nu_hat
    result$vcov <- J %*% result$vcov %*% t(J)
  }

  result
}

.hzr_numeric_grad_weibull <- function(theta, time, status, time_lower = NULL, time_upper = NULL, x = NULL) {
  objective <- function(par) {
    .hzr_logl_weibull(
      theta = par,
      time = time,
      status = status,
      time_lower = time_lower,
      time_upper = time_upper,
      x = x,
      return_gradient = FALSE
    )
  }

  if (requireNamespace("numDeriv", quietly = TRUE)) {
    return(as.numeric(numDeriv::grad(objective, theta)))
  }

  # Fallback central differences if numDeriv is unavailable.
  eps <- 1e-6
  grad <- numeric(length(theta))
  for (i in seq_along(theta)) {
    tp <- theta
    tm <- theta
    tp[i] <- tp[i] + eps
    tm[i] <- tm[i] - eps
    grad[i] <- (objective(tp) - objective(tm)) / (2 * eps)
  }
  grad
}
