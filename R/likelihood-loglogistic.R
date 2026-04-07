#' @importFrom stats optim
#' @keywords internal
NULL

# likelihood-loglogistic.R — Log-logistic parametric hazard likelihood, gradient, and optimizer
#
# MODEL
# -----
# Log-logistic proportional-odds model (commonly used in actuarial work):
#
#   h(t | x)  = α β t^(β−1) / (1 + α t^β exp(η))   hazard
#   H(t | x)  = log(1 + α t^β exp(η))               cumulative hazard
#   S(t | x)  = 1 / (1 + α t^β exp(η))              survival
#
# where α > 0 (scale), β > 0 (shape), η = x β (linear predictor).
#
# When β = 1, the hazard function is monotone decreasing.  When β > 1 it
# has a single interior maximum (unimodal), which makes the log-logistic
# useful for modelling failure rates that rise then fall over time.
#
# THETA LAYOUT
# ------------
#   theta[1]   = log(α)   (unconstrained; α recovered via exp())
#   theta[2]   = log(β)   (unconstrained; β recovered via exp())
#   theta[3:p] = β_coef   (covariate coefficients; unrestricted)
#
# KEY GRADIENT NOTE
# -----------------
# Let prob_weight_i = α t_i^β exp(η_i) / (1 + α t_i^β exp(η_i))
# (probability weight, i.e., the logistic "p" for each observation).
#
#   dL/d(log α) = sum(δ) − sum(prob_weight)                       [n scalar]
#   dL/d(log β) = sum(δ) + sum(δ·log t) − β·sum(log(t)·prob_weight)
#   dL/dβ_j     = −t(X) %*% prob_weight                          [NO δ term]
#
# The β_j gradient has no event-indicator term because the log-logistic
# log-likelihood penalises all observations for the log(1+term) survival
# contribution, not just events.  This differs from, e.g., exponential.
#
# FUNCTIONS
# ---------
#   .hzr_logl_loglogistic()     — log-likelihood (optionally returning gradient)
#   .hzr_gradient_loglogistic() — analytical score vector
#   .hzr_optim_loglogistic()    — unconstrained BFGS wrapper

#' Log-Logistic Parametric Hazard Likelihood
#'
#' Evaluate the log-likelihood and its derivatives for log-logistic hazard models.
#' The log-logistic distribution is more flexible than exponential with separate
#' scale (α) and shape (β) parameters.
#'
#' @keywords internal

#' Log-likelihood for log-logistic hazard with covariates
#'
#' Computes the negative log-likelihood for right-censored data under the log-logistic
#' parametric hazard model with optional linear-predictor covariates.
#'
#' @param theta Vector of parameters:
#'   theta\[1\] = log(α) where α > 0 is the scale parameter
#'   theta\[2\] = log(β) where β > 0 is the shape parameter
#'   theta\[3:length\]: Covariate coefficients (linear on log-scale of cumulative hazard)
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
#' \deqn{\ell(\theta) = \sum_{\delta_i = 1} [\log(\alpha \beta) + (\beta - 1) \log(t_i) - \log(1 + \alpha t_i^{\beta} \exp(\eta_i))]
#'   + \sum_{i} \log(1 + \alpha t_i^{\beta} \exp(\eta_i))}
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
    return_gradient = FALSE) {
  
  n <- length(time)
  
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

  # Exact event: log h + log S
  ll_event <- if (any(idx_event)) {
    sum(
      log_alpha + log_beta + (beta - 1) * log(time[idx_event]) -
        2 * log1p(term_event[idx_event])
    )
  } else {
    0
  }

  # Right-censored: log S = -log(1 + term)
  ll_right <- if (any(idx_right)) {
    -sum(log1p(term_event[idx_right]))
  } else {
    0
  }

  # Left-censored: log F = log(term/(1+term))
  ll_left <- if (any(idx_left)) {
    sum(log_term_u[idx_left] - log1p(term_upper[idx_left]))
  } else {
    0
  }

  # Interval-censored: log(F(u) - F(l)).
  # We compute on probability scale and then log, with positivity guard.
  ll_interval <- if (any(idx_interval)) {
    p_u <- term_upper[idx_interval] / (1 + term_upper[idx_interval])
    p_l <- term_lower[idx_interval] / (1 + term_lower[idx_interval])
    diff_p <- p_u - p_l
    if (any(diff_p <= 0)) return(Inf)
    sum(log(diff_p))
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
      theta, time, status, time_lower, time_upper, x, eta, alpha, beta, log_alpha, log_beta, term_event
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
#'   \eqn{L = \sum(\delta_i * [\log(\alpha) + \log(\beta) + (\beta-1)*\log(t_i)]) - \sum(\log(1 + \alpha*t_i^\beta*\exp(\eta_i)))}
#'
#' Let p_i = α·t_i^β·exp(η_i) / (1 + α·t_i^β·exp(η_i)) = probability of event at t_i
#'
#' Derivatives (using chain rule and reparameterization):
#' dL/d(log α) = sum(δ_i) - α * sum(t_i^β * exp(η_i) / (1 + α*t_i^β*exp(η_i)))
#'             = sum(δ_i) - sum(α * t_i^β * exp(η_i) / (1 + term))
#'
#' dL/d(log β) = sum(δ_i) + sum(δ_i * log(t_i)) - β * sum(α*t_i^β*log(t_i)*exp(η_i)/(1 + α*t_i^β*exp(η_i)))
#'             = sum(δ_i) + sum(δ_i * log(t_i)) - β * sum(t_i^β*log(t_i)/(1 + 1/(α*t_i^β*exp(η_i))))
#'
#' dL/dβ_j = sum(δ_i * x_ij) - sum(α*t_i^β*exp(η_i)*x_ij / (1 + α*t_i^β*exp(η_i)))
#'         = t(X) %*% (δ - p)
#'
#' @noRd
.hzr_gradient_loglogistic <- function(
    theta,
    time,
    status,
  time_lower = NULL,
  time_upper = NULL,
    x = NULL,
    eta = NULL,
    alpha = NULL,
    beta = NULL,
    log_alpha = NULL,
    log_beta = NULL,
    term = NULL) {
  
  # Numerical score currently used for all status patterns to guarantee
  # consistency with the mixed-censoring likelihood expression above.
  .hzr_numeric_grad_loglogistic(theta, time, status, time_lower, time_upper, x)
}

#' Optimize log-logistic hazard likelihood
#'
#' Fits a parametric log-logistic hazard model by maximizing the log-likelihood
#' using L-BFGS-B optimization (unconstrained, via log(alpha) and log(beta) reparameterization).
#'
#' @param time Numeric vector of follow-up times
#' @param status Numeric vector of event indicators (0/1)
#' @param x Optional design matrix of covariates
#' @param theta_start Vector of starting parameter values
#' @param control List of optimization control parameters
#'
#' @return List with:
#'   - \code{par}: Final parameter estimates
#'   - \code{value}: Log-likelihood at solution
#'   - \code{convergence}: Code (0 = success)
#'   - \code{counts}: Number of function/gradient evaluations
#'   - \code{message}: Convergence message
#'   - \code{vcov}: Variance-covariance matrix at solution
#'
#' @noRd
.hzr_optim_loglogistic <- function(
    time,
    status,
  time_lower = NULL,
  time_upper = NULL,
    x = NULL,
    theta_start,
    control = list()) {
  
  # Merge defaults
  control <- utils::modifyList(
    list(
      maxit = 1000,
      reltol = 1e-5,
      abstol = 1e-6,
      method = "BFGS"
    ),
    control
  )
  
  # Set up objective (negative log-likelihood for minimization)
  objective <- function(theta) {
    ll <- -.hzr_logl_loglogistic(
      theta = theta,
      time = time,
      status = status,
      time_lower = time_lower,
      time_upper = time_upper,
      x = x,
      return_gradient = FALSE
    )
    
    if (!is.finite(ll)) return(1e10)
    ll
  }
  
  # Set up gradient
  gradient <- function(theta) {
    # Prevent numerical issues
    if (!all(is.finite(theta))) {
      return(rep(0, length(theta)))
    }
    
    log_alpha <- theta[1]
    log_beta <- theta[2]
    alpha <- exp(log_alpha)
    beta <- exp(log_beta)
    
    if (!is.null(x)) {
      beta_coef <- theta[3:length(theta)]
      eta <- as.numeric(x %*% beta_coef)
    } else {
      eta <- rep(0, length(time))
    }
    
    log_term <- log_alpha + beta * log(time) + eta
    
    # Check for overflow
    if (max(log_term, na.rm = TRUE) > 100) {
      return(rep(0, length(theta)))
    }
    
    term <- exp(log_term)
    
    if (!all(is.finite(term))) {
      return(rep(0, length(theta)))
    }
    
    grad <- tryCatch({
      .hzr_gradient_loglogistic(
        theta = theta,
        time = time,
        status = status,
        time_lower = time_lower,
        time_upper = time_upper,
        x = x,
        eta = eta,
        alpha = alpha,
        beta = beta,
        log_alpha = log_alpha,
        log_beta = log_beta,
        term = term
      )
    }, error = function(e) {
      rep(0, length(theta))
    })
    
    if (!all(is.finite(grad))) {
      grad[!is.finite(grad)] <- 0
    }
    
    # Return negative gradient (for minimization)
    -grad
  }
  
  # Optimize (unconstrained optimization since log(alpha) and log(beta) are unrestricted)
  result <- optim(
    par = theta_start,
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(
      maxit = control$maxit,
      reltol = control$reltol
    ),
    hessian = FALSE
  )
  
  # Post-fit Hessian computation
  hess_result <- tryCatch(
    {
      if (requireNamespace("numDeriv", quietly = TRUE)) {
        numDeriv::hessian(objective, result$par)
      } else {
        NA
      }
    },
    error = function(e) NA
  )
  
  vcov <- if (!anyNA(hess_result)) {
    tryCatch(
      solve(-hess_result),
      error = function(e) {
        warning("Hessian not invertible; standard errors unavailable")
        NA
      }
    )
  } else {
    NA
  }
  
  list(
    par = result$par,
    value = -result$value,  # Restore log-likelihood
    convergence = result$convergence,
    counts = result$counts,
    message = result$message,
    hessian = hess_result,
    vcov = vcov
  )
}

.hzr_numeric_grad_loglogistic <- function(theta, time, status, time_lower = NULL, time_upper = NULL, x = NULL) {
  # Reuse log-likelihood as scalar objective and differentiate numerically.
  objective <- function(par) {
    .hzr_logl_loglogistic(
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
