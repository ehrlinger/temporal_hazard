#' @importFrom stats optim
#' @keywords internal
NULL

# likelihood-exponential.R — Exponential parametric hazard likelihood, gradient, and optimizer
#
# MODEL
# -----
# Constant-hazard (memoryless) proportional-hazards model:
#
#   h(t | x)  = λ exp(η)                    hazard (constant in t)
#   H(t | x)  = λ t exp(η)                  cumulative hazard
#   S(t | x)  = exp(−λ t exp(η))            survival
#
# where λ > 0 (baseline rate), η = x β (linear predictor).
#
# THETA LAYOUT
# ------------
#   theta[1]   = log(λ)   (unconstrained; λ recovered via exp())
#   theta[2:p] = β        (covariate coefficients; unrestricted)
#
# The exponential is a special case of Weibull with ν = 1.  It has no shape
# parameter so theta is one element shorter than Weibull for the same data.
# Unconstrained BFGS is used (no box bounds needed).
#
# FUNCTIONS
# ---------
#   .hzr_logl_exponential()     — log-likelihood (optionally returning gradient)
#   .hzr_gradient_exponential() — analytical score vector
#   .hzr_optim_exponential()    — unconstrained BFGS wrapper

#' Exponential Parametric Hazard Likelihood
#'
#' Evaluate the log-likelihood and its derivatives for exponential hazard models.
#' The exponential is the simplest ALT model: constant hazard λ, no shape parameter.
#'
#' @keywords internal

#' Log-likelihood for exponential hazard with covariates
#'
#' Computes the negative log-likelihood for right-censored data under the exponential
#' parametric hazard model with optional linear-predictor covariates.
#'
#' @param theta Vector of parameters:
#'   theta\[1\] = log(λ) where λ > 0 is the baseline hazard rate
#'   theta\[2:length\]: Covariate coefficients (linear on log-hazard scale)
#'
#' @param time Numeric vector of follow-up times (n)
#' @param status Numeric vector of event indicators: 1 = event, 0 = censored (n)
#' @param time_lower Optional numeric lower bound vector for interval-censored rows.
#'   Defaults to time if NULL.
#' @param time_upper Optional numeric upper bound vector for left/interval-censored rows.
#'   Defaults to time if NULL.
#' @param x Design matrix of covariates (n × p_coef); NULL for no covariates
#' @param return_gradient Logical; if TRUE, attach gradient vector as attribute
#'
#' @return Scalar log-likelihood value. If return_gradient = TRUE, gradient vector
#' is attached as attribute \code{"gradient"}.
#'
#' @details
#' The exponential hazard model is parameterized as:
#'
#' \deqn{h(t | x) = \lambda \exp(\eta)}
#'
#' where:
#' - \eqn{\lambda > 0} is baseline hazard rate
#' - \eqn{\eta = x \beta} is log-relative-hazard (covariates)
#'
#' The survival function is:
#'
#' \deqn{S(t | x) = \exp(-\lambda t \exp(\eta))}
#'
#' The log-likelihood for right-censored data is:
#'
#' \deqn{\ell(\theta) = \sum_{i: \delta_i = 1} [\log\lambda + \eta_i]
#'   - \lambda \sum_i t_i \exp(\eta_i)}
#'
#' Reparameterization: \u03b8\[1\] = log(\u03bb) avoids constrained optimization.
#'
#' Mixed censoring status coding:
#' - 1: exact event at time
#' - 0: right-censored at time
#' - -1: left-censored with upper bound time_upper \(or time\)
#' - 2: interval-censored in the interval \(time_lower, time_upper\)
#'
#' @noRd
.hzr_logl_exponential <- function(
    theta,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    weights = NULL,
    return_gradient = FALSE) {

  n <- length(time)

  # Extract parameters
  log_lambda <- theta[1]
  lambda <- exp(log_lambda)  # Always positive

  # Covariate coefficients (if any)
  if (!is.null(x)) {
    if (is.null(attr(x, "dimnames")[[2]])) {
      colnames(x) <- paste0("beta_", seq_len(ncol(x)))
    }
    beta <- theta[2:length(theta)]
    eta <- as.numeric(x %*% beta)  # n-vector
  } else {
    eta <- rep(0, n)
  }

  # Resolve censoring bounds once so each likelihood term can reference
  # a consistent lower/upper representation regardless of user input style.
  lower <- if (is.null(time_lower)) time else time_lower
  upper <- if (is.null(time_upper)) time else time_upper

  if (any(status %in% c(-1, 2)) && any(upper <= 0)) {
    return(Inf)
  }
  if (any(status == 2) && any(lower[status == 2] >= upper[status == 2])) {
    return(Inf)
  }

  # Cumulative hazards at event/lower/upper times; all terms can be
  # expressed using H(t) for the exponential model.
  cumhaz_event <- lambda * time * exp(eta)
  cumhaz_lower <- lambda * lower * exp(eta)
  cumhaz_upper <- lambda * upper * exp(eta)

  idx_event <- status == 1
  idx_right <- status == 0
  idx_left <- status == -1
  idx_interval <- status == 2

  # Event: log h(t) + log S(t) = [log(lambda)+eta] - H(t)
  ll_event <- if (any(idx_event)) {
    sum((log_lambda + eta[idx_event]) - cumhaz_event[idx_event])
  } else {
    0
  }

  # Right-censored: log S(t) = -H(t)
  ll_right <- if (any(idx_right)) {
    -sum(cumhaz_event[idx_right])
  } else {
    0
  }

  # Left-censored: log F(u) = log(1 - exp(-H(u)))
  ll_left <- if (any(idx_left)) {
    sum(hzr_log1mexp(cumhaz_upper[idx_left]))
  } else {
    0
  }

  # Interval-censored: log(S(l)-S(u))
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

  # If gradient requested, compute score vector
  if (return_gradient) {
    grad <- .hzr_gradient_exponential(
      theta, time, status, time_lower, time_upper, x, eta, cumhaz_event, lambda, log_lambda
    )
    attr(logl, "gradient") <- grad
  }

  logl
}

#' Score vector (gradient of log-likelihood) for exponential
#'
#' Computes the score vector of the exponential log-likelihood w.r.t. all parameters.
#'
#' The log-likelihood is:
#'   \eqn{L = \sum(\delta_i * [\log(\lambda) + \eta_i]) - \sum(\lambda * t_i * \exp(\eta_i))}
#'
#' Derivatives:
#' dL/d(log λ) = sum(δ_i) - sum(λ * t_i * exp(η_i)) = sum(δ_i) - sum(H_i)
#' dL/dβ_j    = sum(δ_i * x_ij) - sum(λ * t_i * exp(η_i) * x_ij) = t(X) %*% (δ - H)
#'
#' @noRd
.hzr_gradient_exponential <- function(
    theta,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    weights = NULL,
    eta = NULL,
    cumhaz = NULL,
    lambda = NULL,
    log_lambda = NULL) {

  n <- length(time)
  p <- length(theta)
  grad <- numeric(p)

  # Numerical fallback for left/interval censoring contributions.
  # Right-censored/event-only rows still use the closed-form score below.
  if (any(status %in% c(-1, 2))) {
    return(.hzr_numeric_grad_exponential(theta, time, status, time_lower, time_upper, x))
  }

  # Recompute if not provided
  if (is.null(eta) || is.null(cumhaz) || is.null(lambda)) {
    log_lambda <- theta[1]
    lambda <- exp(log_lambda)

    if (!is.null(x)) {
      beta <- theta[2:length(theta)]
      eta <- as.numeric(x %*% beta)
    } else {
      eta <- rep(0, n)
    }

    cumhaz <- lambda * time * exp(eta)
  }

  # ===== Gradient w.r.t. log(lambda) =====
  # dL/d(log λ) = sum(δ) - sum(H) = sum(δ) - sum(λ * t * exp(η))
  grad[1] <- sum(status) - sum(cumhaz)

  # ===== Gradient w.r.t. beta (covariate coefficients) =====
  if (p > 1 && !is.null(x)) {
    # dL/dβ = t(X) %*% (δ - H)
    residual <- status - cumhaz
    grad[2:p] <- as.numeric(crossprod(x, residual))
  }

  grad
}

#' @noRd
.hzr_optim_exponential <- function(
    time, status, time_lower = NULL, time_upper = NULL,
    x = NULL, theta_start, weights = NULL, control = list()) {
  .hzr_optim_generic(
    logl_fn = .hzr_logl_exponential,
    gradient_fn = .hzr_gradient_exponential,
    time = time, status = status,
    time_lower = time_lower, time_upper = time_upper,
    x = x, theta_start = theta_start, weights = weights,
    control = control, use_bounds = FALSE
  )
}

.hzr_numeric_grad_exponential <- function(theta, time, status, time_lower = NULL, time_upper = NULL, x = NULL) {
  # Reuse log-likelihood as scalar objective and differentiate numerically.
  objective <- function(par) {
    .hzr_logl_exponential(
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
