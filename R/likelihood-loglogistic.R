#' @importFrom stats optim
#' @keywords internal
NULL

# likelihood-loglogistic.R -- Log-logistic parametric hazard likelihood, gradient, and optimizer
#
# MODEL
# -----
# Log-logistic proportional-odds model (commonly used in actuarial work):
#
#   h(t | x)  = alpha beta t^(beta-1) exp(eta) / (1 + alpha t^beta exp(eta))   hazard
#   H(t | x)  = log(1 + alpha t^beta exp(eta))               cumulative hazard
#   S(t | x)  = 1 / (1 + alpha t^beta exp(eta))              survival
#
# where alpha > 0 (scale), beta > 0 (shape), eta = x beta (linear predictor).
#
# When beta = 1, the hazard function is monotone decreasing.  When beta > 1 it
# has a single interior maximum (unimodal), which makes the log-logistic
# useful for modelling failure rates that rise then fall over time.
#
# THETA LAYOUT
# ------------
#   theta[1]   = log(alpha)   (unconstrained; alpha recovered via exp())
#   theta[2]   = log(beta)   (unconstrained; beta recovered via exp())
#   theta[3:p] = beta_coef   (covariate coefficients; unrestricted)
#
# KEY GRADIENT NOTE
# -----------------
# Let pw_i = term_i / (1 + term_i) where term_i = alpha t_i^beta exp(eta_i)
# and let w_i = (1 + delta_i) * pw_i  (events get weight 2, censored get 1).
#
#   dL/d(log alpha) = sum(delta) - sum(w)
#   dL/d(log beta) = sum(delta) + beta * [sum(delta * log t) - sum(w * log t)]
#   dL/dbeta_j     = t(X) %*% (delta - w)
#
# The (1+delta) weighting arises because events contribute both log h(t) and
# log S(t), each of which depends on log(1 + term).  Censored observations
# contribute only log S(t).
#
# FUNCTIONS
# ---------
#   .hzr_logl_loglogistic()     -- log-likelihood (optionally returning gradient)
#   .hzr_gradient_loglogistic() -- analytical score vector
#   .hzr_optim_loglogistic()    -- unconstrained BFGS wrapper

#' Log-Logistic Parametric Hazard Likelihood
#'
#' Evaluate the log-likelihood and its derivatives for log-logistic hazard models.
#' The log-logistic distribution is more flexible than exponential with separate
#' scale (alpha) and shape (beta) parameters.
#'
#' @keywords internal

#' Log-likelihood for log-logistic hazard with covariates
#'
#' Computes the negative log-likelihood for right-censored data under the log-logistic
#' parametric hazard model with optional linear-predictor covariates.
#'
#' @param theta Vector of parameters:
#'   theta\[1\] = log(alpha) where alpha > 0 is the scale parameter
#'   theta\[2\] = log(beta) where beta > 0 is the shape parameter
#'   theta\[3:length\]: Covariate coefficients (linear on log-scale of cumulative hazard)
#'
#' @param time Numeric vector of follow-up times (n)
#' @param status Numeric vector of event indicators: 1 = event, 0 = censored (n)
#' @param time_lower Optional numeric lower bound vector for interval-censored rows.
#'   Defaults to time if NULL.
#' @param time_upper Optional numeric upper bound vector for left/interval-censored rows.
#'   Defaults to time if NULL.
#' @param x Design matrix of covariates (n x p_coef); NULL for no covariates
#' @param return_gradient Logical; if TRUE, attach gradient vector as attribute
#'
#' @return Scalar log-likelihood value. If return_gradient = TRUE, gradient vector
#' is attached as attribute \code{"gradient"}.
#'
#' @details
#' The log-logistic hazard model is parameterized as:
#'
#' \deqn{h(t | x) = \frac{\alpha \beta t^{\beta - 1}}{1 + \alpha t^{\beta} \exp(\eta)}}
#'
#' where:
#' - \eqn{\alpha > 0} is scale parameter
#' - \eqn{\beta > 0} is shape parameter
#' - \eqn{\eta = x \beta} is covariate effect (linear on log-scale)
#'
#' The survival function is:
#'
#' \deqn{S(t | x) = \frac{1}{1 + \alpha t^{\beta} \exp(\eta)}}
#'
#' The cumulative hazard is:
#'
#' \deqn{H(t | x) = \log(1 + \alpha t^{\beta} \exp(\eta))}
#'
#' The log-likelihood for right-censored data is:
#'
#' \deqn{\ell(\theta) = \sum_{\delta_i = 1} [\log(\alpha \beta) + (\beta - 1)
#'   \log(t_i) - \log(1 + \alpha t_i^{\beta} \exp(\eta_i))] + \sum_{i}
#'   \log(1 + \alpha t_i^{\beta} \exp(\eta_i))}
#'
#' Reparameterization: \u03b8\[1\] = log(\u03b1), \u03b8\[2\] = log(\u03b2) avoids constrained optimization.
#'
#' Mixed censoring status coding:
#' - 1: exact event at time
#' - 0: right-censored at time
#' - -1: left-censored with upper bound time_upper \(or time\)
#' - 2: interval-censored in the interval \(time_lower, time_upper\)
#'
#' @noRd
.hzr_logl_loglogistic <- function(
    theta,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    weights = NULL,
    return_gradient = FALSE) {

  n <- length(time)

  # Default unit weights
  if (is.null(weights)) weights <- rep(1, n)

  # Extract parameters
  log_alpha <- theta[1]
  log_beta <- theta[2]
  alpha <- exp(log_alpha)  # Always positive
  beta <- exp(log_beta)    # Always positive

  # Covariate coefficients (if any)
  if (!is.null(x)) {
    if (is.null(attr(x, "dimnames")[[2]])) {
      colnames(x) <- paste0("beta_", seq_len(ncol(x)))
    }
    beta_coef <- theta[3:length(theta)]
    eta <- as.numeric(x %*% beta_coef)  # n-vector
  } else {
    eta <- rep(0, n)
  }

  lower <- if (is.null(time_lower)) time else time_lower
  upper <- if (is.null(time_upper)) time else time_upper

  if (any(status == 1 & time <= 0)) {
    return(Inf)
  }
  if (any(status %in% c(-1, 2) & upper <= 0)) {
    return(Inf)
  }
  if (any(status == 2) && any(lower[status == 2] >= upper[status == 2])) {
    return(Inf)
  }

  # Build term values for event/lower/upper times. Using log-scale first keeps
  # overflow checks centralized before exponentiation.
  log_term <- log_alpha + beta * log(time) + eta
  log_term_l <- log_alpha + beta * log(lower) + eta
  log_term_u <- log_alpha + beta * log(upper) + eta

  # Check for numerical overflow
  max_log_term <- max(c(log_term, log_term_l, log_term_u), na.rm = TRUE)
  if (max_log_term > 100) {
    return(Inf)
  }

  term_event <- exp(log_term)
  term_lower <- exp(log_term_l)
  term_upper <- exp(log_term_u)

  idx_event <- status == 1
  idx_right <- status == 0
  idx_left <- status == -1
  idx_interval <- status == 2

  # Exact event: w * log f = w * [log h + log S]
  # f(t|x) = alpha beta t^(beta-1) exp(eta) / (1 + alpha t^beta exp(eta))^2
  # log f  = log alpha + log beta + (beta-1) log t + eta - 2 log(1 + term)
  ll_event <- if (any(idx_event)) {
    sum(
      weights[idx_event] *
        (log_alpha + log_beta + (beta - 1) * log(time[idx_event]) +
           eta[idx_event] - 2 * log1p(term_event[idx_event]))
    )
  } else {
    0
  }

  # Right-censored: w * log S = -w * log(1 + term)
  ll_right <- if (any(idx_right)) {
    -sum(weights[idx_right] * log1p(term_event[idx_right]))
  } else {
    0
  }

  # Left-censored: w * log F = w * log(term/(1+term))
  ll_left <- if (any(idx_left)) {
    sum(weights[idx_left] *
          (log_term_u[idx_left] - log1p(term_upper[idx_left])))
  } else {
    0
  }

  # Interval-censored: w * log(F(u) - F(l)).
  # We compute on probability scale and then log, with positivity guard.
  ll_interval <- if (any(idx_interval)) {
    p_u <- term_upper[idx_interval] / (1 + term_upper[idx_interval])
    p_l <- term_lower[idx_interval] / (1 + term_lower[idx_interval])
    diff_p <- p_u - p_l
    if (any(diff_p <= 0)) return(Inf)
    sum(weights[idx_interval] * log(diff_p))
  } else {
    0
  }

  logl <- ll_event + ll_right + ll_left + ll_interval

  if (!is.finite(logl)) {
    return(Inf)
  }

  # If gradient requested, compute score vector
  if (return_gradient) {
    grad <- .hzr_gradient_loglogistic(
      theta = theta, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper, x = x,
      weights = weights,
      eta = eta, alpha = alpha, beta = beta,
      log_alpha = log_alpha, log_beta = log_beta, term = term_event
    )
    attr(logl, "gradient") <- grad
  }

  logl
}

#' Score vector (gradient of log-likelihood) for log-logistic
#'
#' Computes the score vector of the log-logistic log-likelihood w.r.t. all parameters.
#'
#' The log-likelihood is:
#'   \eqn{L = \sum(\delta_i * [\log(\alpha) + \log(\beta) + (\beta-1)
#'   *\log(t_i)]) - \sum(\log(1 + \alpha*t_i^\beta*\exp(\eta_i)))}
#'
#' Let p_i = alpha*t_i^beta*exp(eta_i) / (1 + alpha*t_i^beta*exp(eta_i)) = probability of event at t_i
#'
#' Derivatives (using chain rule and reparameterization):
#' dL/d(log alpha) = sum(delta_i) - alpha * sum(t_i^beta * exp(eta_i) / (1 + alpha*t_i^beta*exp(eta_i)))
#'             = sum(delta_i) - sum(alpha * t_i^beta * exp(eta_i) / (1 + term))
#'
#' dL/d(log beta) = sum(delta_i) + sum(delta_i * log(t_i))
#'                  - beta * sum(alpha*t_i^beta*log(t_i)*exp(eta_i)
#'                                 /(1 + alpha*t_i^beta*exp(eta_i)))
#'                = sum(delta_i) + sum(delta_i * log(t_i))
#'                  - beta * sum(t_i^beta*log(t_i)
#'                                 /(1 + 1/(alpha*t_i^beta*exp(eta_i))))
#'
#' dL/dbeta_j = sum(delta_i * x_ij) - sum(alpha*t_i^beta*exp(eta_i)*x_ij / (1 + alpha*t_i^beta*exp(eta_i)))
#'         = t(X) %*% (delta - p)
#'
#' @noRd
.hzr_gradient_loglogistic <- function(
    theta,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    weights = NULL,
    eta = NULL,
    alpha = NULL,
    beta = NULL,
    log_alpha = NULL,
    log_beta = NULL,
    term = NULL) {

  n <- length(time)
  p <- length(theta)
  grad <- numeric(p)

  # Default unit weights -- keep gradient consistent with the LL when weights
  # are omitted.
  if (is.null(weights)) weights <- rep(1, n)

  # Mixed censoring uses numerical gradient for robustness.
  if (any(status %in% c(-1, 2))) {
    return(.hzr_numeric_grad_loglogistic(theta, time, status, time_lower,
                                          time_upper, x, weights = weights))
  }

  # Recompute if not provided
  if (is.null(term) || is.null(alpha) || is.null(beta)) {
    log_alpha <- theta[1]
    log_beta <- theta[2]
    alpha <- exp(log_alpha)
    beta <- exp(log_beta)

    if (!is.null(x)) {
      beta_coef <- theta[3:length(theta)]
      eta <- as.numeric(x %*% beta_coef)
    } else {
      eta <- rep(0, n)
    }

    term <- alpha * (time ^ beta) * exp(eta)
  }

  # pw_i = term_i / (1 + term_i), wm_i = weights_i * (1 + delta_i) * pw_i
  # (wm is the weighted Mills building block used in every derivative term.)
  pw <- term / (1 + term)
  w_status <- weights * status
  wm <- weights * (1 + status) * pw

  # dL/d(log alpha) = sum(w * delta) - sum(wm)
  grad[1] <- sum(w_status) - sum(wm)

  # dL/d(log beta) = sum(w * delta) + beta * [sum(w * delta * log t) - sum(wm * log t)]
  log_t <- log(time)
  grad[2] <- sum(w_status) +
             beta * (sum(w_status * log_t) - sum(wm * log_t))

  # dL/dbeta_j = t(X) %*% (w * delta - wm)
  if (p > 2 && !is.null(x)) {
    residual <- w_status - wm
    grad[3:p] <- as.numeric(crossprod(x, residual))
  }

  grad
}

#' @noRd
.hzr_optim_loglogistic <- function(
    time, status, time_lower = NULL, time_upper = NULL,
    x = NULL, theta_start, weights = NULL, control = list()) {
  .hzr_optim_generic(
    logl_fn = .hzr_logl_loglogistic,
    gradient_fn = .hzr_gradient_loglogistic,
    time = time, status = status,
    time_lower = time_lower, time_upper = time_upper,
    x = x, theta_start = theta_start, weights = weights,
    control = control, use_bounds = FALSE
  )
}

.hzr_numeric_grad_loglogistic <- function(theta, time, status,
                                            time_lower = NULL, time_upper = NULL,
                                            x = NULL, weights = NULL) {
  # Reuse log-likelihood as scalar objective and differentiate numerically.
  objective <- function(par) {
    .hzr_logl_loglogistic(
      theta = par,
      time = time,
      status = status,
      time_lower = time_lower,
      time_upper = time_upper,
      x = x,
      weights = weights,
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
