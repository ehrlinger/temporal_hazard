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
#'   theta[1:n_shape]: Shape/scale parameters (distribution-specific)
#'   theta[n_shape+1:length]: Covariate coefficients (linear on log-hazard scale)
#'   where n_shape = number of shape parameters (e.g., 2 for Weibull scale/shape)
#'
#' @param time Numeric vector of follow-up times (n)
#' @param status Numeric vector of event indicators: 1 = event, 0 = censored (n)
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
#' @noRd
.hzr_logl_weibull <- function(
    theta,
    time,
    status,
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
  
  # Hazard: h(t) = mu * nu * t^(nu-1) * exp(eta)
  haz <- mu * nu * (time ^ (nu - 1)) * exp(eta)
  
  # Survival (Weibull CDF): S(t) = exp(-(mu*t)^nu * exp(eta))
  cumhaz <- (mu * time) ^ nu * exp(eta)  # cumulative hazard
  surv <- exp(-cumhaz)
  
  # Log-likelihood (negative for minimization)
  logl <- sum(status * log(haz)) - sum(cumhaz)
  
  if (!is.finite(logl)) {
    return(Inf)
  }
  
  # If gradient requested, compute score vector (will be filled in M2-full port)
  if (return_gradient) {
    grad <- .hzr_gradient_weibull(
      theta, time, status, x, eta, cumhaz, haz, mu, nu, n_shape
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
#' @noRd
.hzr_gradient_weibull <- function(
    theta,
    time,
    status,
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
  
  # Placeholder: full implementation deferred to M2 likelihood porting
  # This is a stub showing the structure. Real implementation requires
  # symbolic or numerical differentiation of .hzr_logl_weibull
  
  grad
}

#' Optimize Weibull hazard likelihood
#'
#' Fits a parametric Weibull hazard model by maximizing the log-likelihood
#' using trust-region or line-search optimization.
#'
#' @param time Numeric vector of follow-up times
#' @param status Numeric vector of event indicators (0/1)
#' @param x Optional design matrix of covariates
#' @param theta_start Vector of starting parameter values
#' @param control List of optimization control parameters:
#'   - \code{maxit}: Maximum iterations (default 1000)
#'   - \code{reltol}: Relative tolerance for convergence (default 1e-5)
#'   - \code{abstol}: Absolute tolerance for gradient norm (default 1e-6)
#'   - \code{method}: "bfgs" or "nm" (Nelder-Mead; default "bfgs")
#'   - \code{condition}: Condition number control (default 14)
#'
#' @return List with:
#'   - \code{par}: Final parameter estimates
#'   - \code{value}: Log-likelihood at solution
#'   - \code{convergence}: Code (0 = success)
#'   - \code{counts}: Number of function/gradient evaluations
#'   - \code{message}: Convergence message
#'   - \code{hessian}: Hessian matrix at solution (for SE computation)
#'
#' @noRd
.hzr_optim_weibull <- function(
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
      method = "BFGS",
      condition = 14
    ),
    control
  )
  
  # Normalize method name
  method <- toupper(control$method)
  if (!method %in% c("NELDER-MEAD", "BFGS", "CG", "L-BFGS-B", "SANN", "BRENT")) {
    stop("method must be one of: BFGS, Nelder-Mead, CG, L-BFGS-B, SANN, Brent", call. = FALSE)
  }
  
  # Set up objective (negative log-likelihood for minimization)
  objective <- function(theta) {
    # Safeguard: parameters must be valid
    if (any(theta[1:2] <= 0)) return(1e10)  # mu, nu > 0
    
    ll <- -.hzr_logl_weibull(
      theta = theta,
      time = time,
      status = status,
      x = x,
      return_gradient = FALSE
    )
    
    # Return large penalty if likelihood is not finite
    if (!is.finite(ll)) return(1e10)
    ll
  }
  
  # Optimize using stats::optim with bounds checking via L-BFGS-B
  # (allows box constraints on parameters: mu > 0, nu > 0)
  lower_bounds <- rep(1e-6, length(theta_start))  # All params > 1e-6
  upper_bounds <- rep(Inf, length(theta_start))
  
  result <- optim(
    par = theta_start,
    fn = objective,
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds,
    control = list(
      maxit = control$maxit,
      factr = 1 / control$reltol,  # L-BFGS-B uses factr instead of reltol
      pgtol = control$abstol       # L-BFGS-B uses pgtol instead of abstol
    ),
    hessian = FALSE  # Disable numerical Hessian; we'll compute it separately if needed
  )
  
  # Post-fit Hessian computation (more stable than from optim)
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
    value = -result$value,  # Restore log-likelihood (un-negate)
    convergence = result$convergence,
    counts = result$counts,
    message = result$message,
    hessian = hess_result,
    vcov = vcov
  )
}
