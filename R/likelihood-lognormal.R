#' @importFrom stats dnorm pnorm optim
#' @keywords internal
NULL

# likelihood-lognormal.R — Log-normal parametric hazard likelihood, gradient, and optimizer
#
# MODEL
# -----
# Accelerated Failure Time (AFT) log-normal model:
#
#   log(T_i) ~ Normal(η_i, σ²)   where η_i = μ + x_i β
#
#   z_i       = (log(t_i) − η_i) / σ          standardised residual
#   S(t | x)  = Φ(−z_i)                        survival  (normal CDF)
#   f(t | x)  = φ(z_i) / (σ t_i)              density
#   h(t | x)  = φ(z_i) / (σ t_i Φ(−z_i))     hazard
#   H(t | x)  = −log Φ(−z_i)                  cumulative hazard
#
# where μ ∈ ℝ (location), σ > 0 (scale), η_i = μ + x_i β (AFT predictor).
#
# AFT vs PH DISTINCTION
# ---------------------
# Unlike Weibull / exponential / log-logistic (PH models), covariates in the
# AFT parameterisation *shift the log-time distribution*.  The linear predictor
# adds to μ (location), not to log-hazard.  Consequently:
#
#   • predict() "linear_predictor" returns β (same interface as PH)
#   • predict() "survival" returns Φ(−z) directly (NOT exp(−H))
#   • predict() "cumulative_hazard" returns −log Φ(−z)
#
# See the EXCEPTION note in hazard_api.R's predict() section banner.
#
# THETA LAYOUT
# ------------
#   theta[1]   = μ         (location; unrestricted)
#   theta[2]   = log(σ)   (unconstrained; σ recovered via exp())
#   theta[3:p] = β        (covariate AFT coefficients; unrestricted)
#
# GRADIENT (let w_i = φ(z_i)/Φ(−z_i), the inverse Mills ratio)
# ---------------------------------------------------------------
#   dL/dμ          = (1/σ) · [sum(δ z) + sum((1−δ) w)]
#   dL/d(log σ)    = sum(δ(z²−1)) + sum((1−δ) w z)
#   dL/dβ_j        = (1/σ) · sum([δ z + (1−δ) w] · x_j)
#
# FUNCTIONS
# ---------
#   .hzr_logl_lognormal()     — log-likelihood (optionally returning gradient)
#   .hzr_gradient_lognormal() — analytical score vector
#   .hzr_optim_lognormal()    — unconstrained BFGS wrapper

#' Log-Normal Parametric Hazard Likelihood
#'
#' Evaluate the log-likelihood and its derivatives for log-normal hazard models.
#' The log-normal is an AFT (accelerated failure time) model with separate
#' location (μ) and scale (σ) parameters.
#'
#' @keywords internal

#' Log-likelihood for log-normal hazard with covariates
#'
#' Computes the log-likelihood for right-censored data under the log-normal
#' parametric hazard model with optional linear-predictor covariates.
#'
#' @param theta Vector of parameters:
#'   theta\[1\] = μ (location, unrestricted)
#'   theta\[2\] = log(σ) where σ > 0 is the scale parameter
#'   theta\[3:length\]: Covariate coefficients (AFT parameterization: add to location)
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
#' The log-normal hazard model uses AFT parameterization:
#'
#' z_i = (log(t_i) - η_i) / σ, where η_i = μ + x_i β
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

  # Extract parameters
  mu <- theta[1]          # Location (unrestricted)
  log_sigma <- theta[2]   # Log-scale (unrestricted)
  sigma <- exp(log_sigma) # Always positive

  # Covariate coefficients (if any) — AFT parameterization: add to location
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

  # Event: log f(t)
  ll_event <- if (any(idx_event)) {
    sum(log_phi_z[idx_event] - log_sigma - log_t[idx_event])
  } else {
    0
  }

  # Right-censored: log S(t)
  ll_right <- if (any(idx_right)) {
    sum(log_surv[idx_right])
  } else {
    0
  }

  # Left-censored: log F(u)
  ll_left <- if (any(idx_left)) {
    z_u <- (log(upper[idx_left]) - eta[idx_left]) / sigma
    sum(pnorm(z_u, log.p = TRUE))
  } else {
    0
  }

  # Interval-censored: log(F(u) - F(l)) with positivity guard.
  ll_interval <- if (any(idx_interval)) {
    z_l <- (log(lower[idx_interval]) - eta[idx_interval]) / sigma
    z_u <- (log(upper[idx_interval]) - eta[idx_interval]) / sigma
    f_l <- pnorm(z_l)
    f_u <- pnorm(z_u)
    diff_f <- f_u - f_l
    if (any(diff_f <= 0)) return(Inf)
    sum(log(diff_f))
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
#' Let z_i = (log(t_i) - η_i) / σ, where η_i = μ + x_i β
#' Let w_i = φ(z_i) / Φ(-z_i) (inverse Mills ratio)
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

  # Mixed-censoring score uses numerical differentiation for robustness.
  if (any(status %in% c(-1, 2))) {
    return(.hzr_numeric_grad_lognormal(theta, time, status, time_lower, time_upper, x))
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

  # Inverse Mills ratio: w_i = φ(z_i) / Φ(-z_i)
  # Computed in log scale for numerical stability, then exponentiating
  log_w <- log_phi_z - log_surv
  w <- exp(log_w)  # Safe: log_surv can be -Inf when Φ(-z) → 0 (z → +∞)
  # Handle case where log_surv = -Inf (z very large → all censored obs have w → ∞)
  # But this only matters for censored obs. For events, w is not used in dL/dμ.
  # Clamp to prevent Inf * 0 issues.
  w <- pmin(w, 1e6)

  censored <- 1 - status  # (1 - δ_i)

  # ===== Gradient w.r.t. μ =====
  # dL/dμ = (1/σ) * [sum(δ_i * z_i) + sum((1-δ_i) * w_i)]
  grad[1] <- (sum(status * z) + sum(censored * w)) / sigma

  # ===== Gradient w.r.t. log(σ) =====
  # dL/d(log σ) = sum(δ_i * (z_i² - 1)) + sum((1-δ_i) * w_i * z_i)
  grad[2] <- sum(status * (z^2 - 1)) + sum(censored * w * z)

  # ===== Gradient w.r.t. covariate β =====
  # dL/dβ_j = (1/σ) * sum([δ_i * z_i + (1-δ_i) * w_i] * x_ij)
  if (p > 2 && !is.null(x)) {
    score_i <- (status * z + censored * w) / sigma
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

.hzr_numeric_grad_lognormal <- function(theta, time, status, time_lower = NULL, time_upper = NULL, x = NULL) {
  # Reuse log-likelihood as scalar objective and differentiate numerically.
  objective <- function(par) {
    .hzr_logl_lognormal(
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
