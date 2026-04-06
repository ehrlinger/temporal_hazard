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
#'   theta[1] = log(λ) where λ > 0 is the baseline hazard rate
#'   theta[2:length]: Covariate coefficients (linear on log-hazard scale)
#'
#' @param time Numeric vector of follow-up times (n)
#' @param status Numeric vector of event indicators: 1 = event, 0 = censored (n)
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
#' Reparameterization: θ[1] = log(λ) avoids constrained optimization.
#'
#' @noRd
.hzr_logl_exponential <- function(
    theta,
    time,
    status,
    x = NULL,
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
  
  # Hazard: h(t) = lambda * exp(eta)
  # Survival: S(t) = exp(-lambda * t * exp(eta))
  # Cumulative hazard: H(t) = lambda * t * exp(eta)
  
  cumhaz <- lambda * time * exp(eta)  # cumulative hazard
  
  # Log-likelihood:
  # ℓ = sum(δ_i * [log(λ) + η_i]) - sum(λ * t_i * exp(η_i))
  logl <- sum(status * (log_lambda + eta)) - sum(cumhaz)
  
  if (!is.finite(logl)) {
    return(Inf)
  }
  
  # If gradient requested, compute score vector
  if (return_gradient) {
    grad <- .hzr_gradient_exponential(
      theta, time, status, x, eta, cumhaz, lambda, log_lambda
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
#'   L = sum(δ_i * [log(λ) + η_i]) - sum(λ * t_i * exp(η_i))
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
    x = NULL,
    eta = NULL,
    cumhaz = NULL,
    lambda = NULL,
    log_lambda = NULL) {
  
  n <- length(time)
  p <- length(theta)
  grad <- numeric(p)
  
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

#' Optimize exponential hazard likelihood
#'
#' Fits a parametric exponential hazard model by maximizing the log-likelihood
#' using L-BFGS-B optimization (unconstrained, via log(lambda) reparameterization).
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
.hzr_optim_exponential <- function(
    time,
    status,
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
    ll <- -.hzr_logl_exponential(
      theta = theta,
      time = time,
      status = status,
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
    
    log_lambda <- theta[1]
    lambda <- exp(log_lambda)
    
    if (!is.null(x)) {
      beta <- theta[2:length(theta)]
      eta <- as.numeric(x %*% beta)
    } else {
      eta <- rep(0, length(time))
    }
    
    cumhaz <- lambda * time * exp(eta)
    
    if (!all(is.finite(cumhaz))) {
      return(rep(0, length(theta)))
    }
    
    grad <- tryCatch({
      .hzr_gradient_exponential(
        theta = theta,
        time = time,
        status = status,
        x = x,
        eta = eta,
        cumhaz = cumhaz,
        lambda = lambda,
        log_lambda = log_lambda
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
  
  # Optimize (unconstrained optimization since log(lambda) is unrestricted)
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
