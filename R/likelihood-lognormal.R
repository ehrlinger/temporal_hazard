#' @importFrom stats dnorm pnorm optim
#' @keywords internal
NULL

# likelihood-lognormal.R -- Log-normal parametric hazard likelihood, gradient, and optimizer
#
# MODEL
# -----
# Accelerated Failure Time (AFT) log-normal model:
#
#   log(T_i) ~ Normal(eta_i, sigma^2)   where eta_i = mu + x_i beta
#
#   z_i       = (log(t_i) - eta_i) / sigma          standardised residual
#   S(t | x)  = Phi(-z_i)                        survival  (normal CDF)
#   f(t | x)  = phi(z_i) / (sigma t_i)              density
#   h(t | x)  = phi(z_i) / (sigma t_i Phi(-z_i))     hazard
#   H(t | x)  = -log Phi(-z_i)                  cumulative hazard
#
# where mu in R (location), sigma > 0 (scale), eta_i = mu + x_i beta (AFT predictor).
#
# AFT vs PH DISTINCTION
# ---------------------
# Unlike Weibull / exponential / log-logistic (PH models), covariates in the
# AFT parameterisation *shift the log-time distribution*.  The linear predictor
# adds to mu (location), not to log-hazard.  Consequently:
#
#   * predict() "linear_predictor" returns beta (same interface as PH)
#   * predict() "survival" returns Phi(-z) directly (NOT exp(-H))
#   * predict() "cumulative_hazard" returns -log Phi(-z)
#
# See the EXCEPTION note in hazard_api.R's predict() section banner.
#
# THETA LAYOUT
# ------------
#   theta[1]   = mu         (location; unrestricted)
#   theta[2]   = log(sigma)   (unconstrained; sigma recovered via exp())
#   theta[3:p] = beta        (covariate AFT coefficients; unrestricted)
#
# GRADIENT (let w_i = phi(z_i)/Phi(-z_i), the inverse Mills ratio)
# ---------------------------------------------------------------
#   dL/dmu          = (1/sigma) * [sum(delta z) + sum((1-delta) w)]
#   dL/d(log sigma)    = sum(delta(z^2-1)) + sum((1-delta) w z)
#   dL/dbeta_j        = (1/sigma) * sum([delta z + (1-delta) w] * x_j)
#
# FUNCTIONS
# ---------
#   .hzr_logl_lognormal()     -- log-likelihood (optionally returning gradient)
#   .hzr_gradient_lognormal() -- analytical score vector
#   .hzr_optim_lognormal()    -- unconstrained BFGS wrapper

#' Log-Normal Parametric Hazard Likelihood
#'
#' Evaluate the log-likelihood and its derivatives for log-normal hazard models.
#' The log-normal is an AFT (accelerated failure time) model with separate
#' location (mu) and scale (sigma) parameters.
#'
#' @keywords internal

#' Log-likelihood for log-normal hazard with covariates
#'
#' Computes the log-likelihood for right-censored data under the log-normal
#' parametric hazard model with optional linear-predictor covariates.
#'
#' @param theta Vector of parameters:
#'   theta\[1\] = mu (location, unrestricted)
#'   theta\[2\] = log(sigma) where sigma > 0 is the scale parameter
#'   theta\[3:length\]: Covariate coefficients (AFT parameterization: add to location)
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
#' The log-normal hazard model uses AFT parameterization:
#'
#' z_i = (log(t_i) - eta_i) / sigma, where eta_i = mu + x_i beta
#'
#' Survival function:
#' \deqn{S(t | x) = \Phi(-z_i)}
#'
#' Density:
#' \deqn{f(t | x) = \phi(z_i) / (\sigma t_i)}
#'
#' Hazard:
#' \deqn{h(t | x) = \phi(z_i) / (\sigma t_i \Phi(-z_i))}
#'
#' The log-likelihood for right-censored data is:
#' \deqn{\ell(\theta) = \sum_{\delta_i=1} [\log \phi(z_i) - \log \sigma - \log t_i]
#'   + \sum_{\delta_i=0} \log \Phi(-z_i)}
#'
#' Reparameterization: \u03b8\[1\] = \u03bc, \u03b8\[2\] = log(\u03c3) avoids constraints.
#'
#' Mixed censoring status coding:
#' - 1: exact event at time
#' - 0: right-censored at time
#' - -1: left-censored with upper bound time_upper \(or time\)
#' - 2: interval-censored in the interval \(time_lower, time_upper\)
#'
#' @noRd
.hzr_logl_lognormal <- function(
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
  mu <- theta[1]          # Location (unrestricted)
  log_sigma <- theta[2]   # Log-scale (unrestricted)
  sigma <- exp(log_sigma) # Always positive

  # Covariate coefficients (if any) -- AFT parameterization: add to location
  if (!is.null(x)) {
    if (is.null(attr(x, "dimnames")[[2]])) {
      colnames(x) <- paste0("beta_", seq_len(ncol(x)))
    }
    beta_coef <- theta[3:length(theta)]
    eta <- mu + as.numeric(x %*% beta_coef)  # Total linear predictor (location)
  } else {
    eta <- rep(mu, n)
  }

  # Normalize censoring bounds to lower/upper vectors for unified formulas.
  lower <- if (is.null(time_lower)) time else time_lower
  upper <- if (is.null(time_upper)) time else time_upper

  if (any(status %in% c(1, 0) & time <= 0)) return(Inf)
  if (any(status == -1 & upper <= 0)) return(Inf)
  if (any(status == 2 & (lower <= 0 | upper <= 0))) return(Inf)
  if (any(status == 2) && any(lower[status == 2] >= upper[status == 2])) return(Inf)

  # Event-time terms are based on z at observed `time`.
  log_t <- log(time)
  z <- (log_t - eta) / sigma
  log_phi_z <- dnorm(z, log = TRUE)
  log_surv <- pnorm(-z, log.p = TRUE)

  idx_event <- status == 1
  idx_right <- status == 0
  idx_left <- status == -1
  idx_interval <- status == 2

  # Event: w * log f(t)
  ll_event <- if (any(idx_event)) {
    sum(weights[idx_event] *
          (log_phi_z[idx_event] - log_sigma - log_t[idx_event]))
  } else {
    0
  }

  # Right-censored: w * log S(t)
  ll_right <- if (any(idx_right)) {
    sum(weights[idx_right] * log_surv[idx_right])
  } else {
    0
  }

  # Left-censored: w * log F(u)
  ll_left <- if (any(idx_left)) {
    z_u <- (log(upper[idx_left]) - eta[idx_left]) / sigma
    sum(weights[idx_left] * pnorm(z_u, log.p = TRUE))
  } else {
    0
  }

  # Interval-censored: w * log(F(u) - F(l)) with positivity guard.
  ll_interval <- if (any(idx_interval)) {
    z_l <- (log(lower[idx_interval]) - eta[idx_interval]) / sigma
    z_u <- (log(upper[idx_interval]) - eta[idx_interval]) / sigma
    f_l <- pnorm(z_l)
    f_u <- pnorm(z_u)
    diff_f <- f_u - f_l
    if (any(diff_f <= 0)) return(Inf)
    sum(weights[idx_interval] * log(diff_f))
  } else {
    0
  }

  logl <- ll_event + ll_right + ll_left + ll_interval

  if (!is.finite(logl)) {
    return(Inf)
  }

  if (return_gradient) {
    grad <- .hzr_gradient_lognormal(
      theta = theta, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper, x = x,
      weights = weights,
      eta = eta, sigma = sigma, log_sigma = log_sigma,
      z = z, log_phi_z = log_phi_z, log_surv = log_surv
    )
    attr(logl, "gradient") <- grad
  }

  logl
}

#' Score vector (gradient of log-likelihood) for log-normal
#'
#' Computes the score vector of the log-normal log-likelihood w.r.t. all parameters.
#'
#' Let z_i = (log(t_i) - eta_i) / sigma, where eta_i = mu + x_i beta
#' Let w_i = phi(z_i) / Phi(-z_i) (inverse Mills ratio)
#'
#' Derivatives:
#' \eqn{dL/d\mu = (1/\sigma) * [\sum(\delta_i * z_i) + \sum((1 - \delta_i) * w_i)]}
#' \eqn{dL/d(\log \sigma) = \sum(\delta_i * (z_i^2 - 1)) + \sum((1 - \delta_i) * w_i * z_i)}
#' \eqn{dL/d\beta_j = (1/\sigma) * \sum([\delta_i * z_i + (1 - \delta_i) * w_i] * x_{ij})}
#'
#' @noRd
.hzr_gradient_lognormal <- function(
    theta,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    weights = NULL,
    eta = NULL,
    sigma = NULL,
    log_sigma = NULL,
    z = NULL,
    log_phi_z = NULL,
    log_surv = NULL) {

  n <- length(time)
  p <- length(theta)
  grad <- numeric(p)

  # Default unit weights -- keep gradient consistent with the LL when weights
  # are omitted.
  if (is.null(weights)) weights <- rep(1, n)

  # Mixed-censoring score uses numerical differentiation for robustness.
  if (any(status %in% c(-1, 2))) {
    return(.hzr_numeric_grad_lognormal(theta, time, status, time_lower,
                                        time_upper, x, weights = weights))
  }

  # Recompute if not provided
  if (is.null(z) || is.null(sigma)) {
    mu <- theta[1]
    log_sigma <- theta[2]
    sigma <- exp(log_sigma)

    if (!is.null(x)) {
      beta_coef <- theta[3:length(theta)]
      eta <- mu + as.numeric(x %*% beta_coef)
    } else {
      eta <- rep(mu, n)
    }

    log_t <- log(time)
    z <- (log_t - eta) / sigma
    log_phi_z <- dnorm(z, log = TRUE)
    log_surv <- pnorm(-z, log.p = TRUE)
  }

  # Inverse Mills ratio: mills_i = phi(z_i) / Phi(-z_i)
  # Computed in log scale for numerical stability, then exponentiating.
  log_mills <- log_phi_z - log_surv
  mills <- exp(log_mills)  # Safe: log_surv can be -Inf when Phi(-z) -> 0 (z -> +Inf)
  # Clamp to prevent Inf * 0 issues when z is very large.
  mills <- pmin(mills, 1e6)

  censored <- 1 - status  # (1 - delta_i)

  # Weighted per-row building blocks. Every term below is the unweighted form
  # with `delta`, `(1-delta)*mills` etc. multiplied by `weights`.
  w_status <- weights * status
  w_censmills <- weights * censored * mills

  # ===== Gradient w.r.t. mu =====
  # dL/dmu = (1/sigma) * sum(w * [delta * z + (1-delta) * mills])
  grad[1] <- (sum(w_status * z) + sum(w_censmills)) / sigma

  # ===== Gradient w.r.t. log(sigma) =====
  # dL/d(log sigma) = sum(w * delta * (z^2 - 1)) + sum(w * (1-delta) * mills * z)
  grad[2] <- sum(w_status * (z^2 - 1)) + sum(w_censmills * z)

  # ===== Gradient w.r.t. covariate beta =====
  # dL/dbeta_j = (1/sigma) * sum(w * [delta * z + (1-delta) * mills] * x_ij)
  if (p > 2 && !is.null(x)) {
    score_i <- (w_status * z + w_censmills) / sigma
    grad[3:p] <- as.numeric(crossprod(x, score_i))
  }

  grad
}

#' @noRd
.hzr_optim_lognormal <- function(
    time, status, time_lower = NULL, time_upper = NULL,
    x = NULL, theta_start, weights = NULL, control = list()) {
  .hzr_optim_generic(
    logl_fn = .hzr_logl_lognormal,
    gradient_fn = .hzr_gradient_lognormal,
    time = time, status = status,
    time_lower = time_lower, time_upper = time_upper,
    x = x, theta_start = theta_start, weights = weights,
    control = control, use_bounds = FALSE
  )
}

.hzr_numeric_grad_lognormal <- function(theta, time, status,
                                          time_lower = NULL, time_upper = NULL,
                                          x = NULL, weights = NULL) {
  # Reuse log-likelihood as scalar objective and differentiate numerically.
  objective <- function(par) {
    .hzr_logl_lognormal(
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
