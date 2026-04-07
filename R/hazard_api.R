#' @importFrom stats optim pnorm
#' @importFrom survival Surv
#' @keywords internal
NULL

# hazard_api.R — Primary user-facing API for TemporalHazard
#
# This file contains the main entry points:
#   hazard()         — construct and optionally fit a parametric hazard model
#   predict.hazard() — generate predictions from a fitted hazard object
#   print.hazard()   — compact S3 print method
#   coef.hazard()    — extract fitted parameter vector
#   vcov.hazard()    — extract variance-covariance matrix
#
# DISTRIBUTION PARAMETERIZATION CONVENTIONS
# -----------------------------------------
# Each distribution stores a flattened theta vector.  Shape/scale parameters
# always come first; covariate coefficients (β) follow:
#
#   Weibull:     theta = [mu, nu, β₁, β₂, ...]      mu>0, nu>0
#   Exponential: theta = [log(λ), β₁, β₂, ...]      λ>0 via exp()
#   Loglogistic: theta = [log(α), log(β), β₁, ...]   α>0, β>0 via exp()
#   Lognormal:   theta = [μ, log(σ), β₁, β₂, ...]   σ>0 via exp(); AFT model
#
# ADDING A NEW DISTRIBUTION
# -------------------------
# 1. Create R/likelihood-<dist>.R with .hzr_logl_<dist>(),
#    .hzr_gradient_<dist>(), and .hzr_optim_<dist>().
# 2. Add an else-if branch in hazard() below under "Distribution dispatch".
# 3. Add an else-if branch in predict.hazard() under "Dispatch by distribution".
# 4. Update the supported-distributions guard in predict.hazard().
# 5. Update the coefficient-extraction block in predict.hazard() for
#    linear_predictor / hazard prediction types.
# 6. Write tests in tests/testthat/test-<dist>-dist.R.
# 7. Generate a golden fixture via .hzr_create_<dist>_golden_fixture().

#' Build and optionally fit a hazard model
#'
#' Creates a `hazard` object and optionally fits it via maximum likelihood.
#' This mirrors the argument-oriented workflow of the legacy HAZARD C/SAS
#' implementation: supply starting values in `theta` and the function will
#' optimize to produce fitted estimates.
#'
#' @param time Numeric follow-up time vector.
#' @param status Numeric or logical event indicator vector.
#' @param time_lower Optional numeric lower bound vector for censoring intervals.
#'   Used when `status == 2` (interval-censored); defaults to `time` if NULL.
#' @param time_upper Optional numeric upper bound vector for censoring intervals.
#'   Used when `status %in% c(-1, 2)`; defaults to `time` if NULL.
#' @param x Optional design matrix (or data frame coercible to matrix).
#' @param formula Optional formula of the form `Surv(time, status) ~ predictors`.
#'   When provided, overrides direct time/status/x arguments and extracts from data.
#'   Example: `hazard(Surv(time, status) ~ x1 + x2, data = df, dist = "weibull", fit = TRUE)`.
#' @param data Optional data frame containing variables referenced in formula.
#' @param time_windows Optional numeric vector of strictly positive cut points for
#'   piecewise time-varying coefficients. When provided, each predictor column in
#'   `x` is expanded into one column per time window so each window gets its own
#'   coefficient.
#' @param theta Optional numeric coefficient vector (starting values for optimization).
#' @param dist Character baseline distribution label (default "weibull").
#' @param fit Logical; if TRUE and theta is provided, fit the model via ML (default TRUE).
#' @param control Named list of control options (see Details).
#' @param ... Additional named arguments retained for parity with legacy calling
#'   conventions.
#'
#' @details
#' Control parameters:
#' - `maxit`: Maximum iterations (default 1000)
#' - `reltol`: Relative parameter change tolerance (default 1e-5)
#' - `abstol`: Absolute gradient norm tolerance (default 1e-6)
#' - `method`: Optimization method: "bfgs" or "nm" (default "bfgs")
#' - `condition`: Condition number control (default 14)
#' - `nocov`, `nocor`: Suppress covariance/correlation output (legacy; no-op in M2)
#'
#' Censoring status coding:
#' - 1: Exact event at time
#' - 0: Right-censored at time
#' - -1: Left-censored with upper bound at time_upper \(or time\)
#' - 2: Interval-censored in the interval \(time_lower, time_upper\)
#'
#' Time-varying coefficients:
#' - If `time_windows` is supplied, predictors are expanded to piecewise window
#'   interactions so each window has its own coefficient vector.
#' - This is implemented as design-matrix expansion, so the existing likelihood
#'   engines remain unchanged.
#'
 #' @export
hazard <- function(formula = NULL,
                   data = NULL,
                   time = NULL,
                   status = NULL,
                   time_lower = NULL,
                   time_upper = NULL,
                   x = NULL,
                   time_windows = NULL,
                   theta = NULL,
                   dist = "weibull",
                   fit = FALSE,
                   control = list(),
                   ...) {
  # Formula dispatch: if formula is provided, parse it and extract time/status/x from data
  if (!is.null(formula)) {
    if (is.null(data)) {
      stop("'data' is required when 'formula' is provided.", call. = FALSE)
    }
    parsed <- .hzr_parse_formula(formula = formula, data = data)
    time <- parsed$time
    status <- parsed$status
    time_lower <- parsed$time_lower
    time_upper <- parsed$time_upper
    x <- parsed$x
  }

  # After formula dispatch, require time and status
  if (is.null(time) || is.null(status)) {
    stop("'time' and 'status' are required (either directly or via 'formula').", call. = FALSE)
  }
  if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0)) {
    stop("'time' must be a numeric vector of finite non-negative values.", call. = FALSE)
  }

  n <- length(time)
  if (length(status) != n) {
    stop("'status' must have the same length as 'time'.", call. = FALSE)
  }

  # Convert Surv object status to numeric if needed (after formula parsing)
  if (inherits(status, "Surv")) {
    status <- unclass(status)[, 2L]
  }

  # Optional censoring bounds:
  # - status = 1 (event): uses observed event time in `time`
  # - status = 0 (right-censored): censoring time in `time` (or `time_lower`)
  # - status = -1 (left-censored): upper bound in `time` (or `time_upper`)
  # - status = 2 (interval-censored): [time_lower, time_upper] required
  if (!is.null(time_lower)) {
    if (!is.numeric(time_lower) || length(time_lower) != n || any(!is.finite(time_lower)) || any(time_lower < 0)) {
      stop("'time_lower' must be a numeric vector of finite non-negative values matching length(time).", call. = FALSE)
    }
  }

  if (!is.null(time_upper)) {
    if (!is.numeric(time_upper) || length(time_upper) != n || any(!is.finite(time_upper)) || any(time_upper < 0)) {
      stop("'time_upper' must be a numeric vector of finite non-negative values matching length(time).", call. = FALSE)
    }
  }

  if (!is.null(x)) {
    x <- .hzr_as_design_matrix(x, n = n)
  }

  if (!is.null(time_windows)) {
    # We use cut points to define window-specific coefficient blocks.
    if (!is.numeric(time_windows) || any(!is.finite(time_windows)) || anyDuplicated(time_windows)) {
      stop("'time_windows' must be a numeric vector of unique finite cut points.", call. = FALSE)
    }
    time_windows <- sort(time_windows)
    if (any(time_windows <= 0)) {
      stop("'time_windows' cut points must be strictly positive.", call. = FALSE)
    }
    if (is.null(x)) {
      stop("'time_windows' requires predictor matrix 'x'.", call. = FALSE)
    }
  }

  x_fit <- x
  if (!is.null(time_windows) && !is.null(x)) {
    # Expand X -> [X_w1 | X_w2 | ...] where each row is active only in its window.
    x_fit <- .hzr_expand_time_varying_design(x = x, time = time, time_windows = time_windows)
  }

  if (!is.null(theta)) {
    if (!is.numeric(theta) || any(!is.finite(theta))) {
      stop("'theta' must be a finite numeric vector when provided.", call. = FALSE)
    }

    # If x exists, theta must include coefficients for all variates
    # theta = [shape parms ... | covariate coefficients ...]
    # For now, assume theta length determines whether we expect x
    required_coef <- if (is.null(x_fit)) 0L else ncol(x_fit)
    if (!is.null(x_fit) && length(theta) < required_coef) {
      stop("'theta' length must be >= number of required coefficients (", required_coef, ").", call. = FALSE)
    }
  }

  if (!is.character(dist) || length(dist) != 1 || !nzchar(dist)) {
    stop("'dist' must be a non-empty character scalar.", call. = FALSE)
  }

  if (!is.list(control)) {
    stop("'control' must be a list.", call. = FALSE)
  }

  if (any(status %in% c(-1, 2))) {
    if (is.null(time_upper)) {
      warning("'time_upper' not provided; using 'time' as upper bound for left/interval-censored rows.")
    }
    if (any(status == 2) && is.null(time_lower)) {
      warning("'time_lower' not provided; using 'time' as lower bound for interval-censored rows.")
    }
  }

  # fit_state holds the result of optimization (or just starting values if fit=FALSE).
  # Fields:
  #   theta     — parameter vector (starting values before fit; MLE estimates after)
  #   converged — TRUE if optimizer reported convergence (code 0); NA if not fitted
  #   objective — log-likelihood at theta; NA_real_ if not fitted
  #   se        — standard errors (sqrt of vcov diagonal); NULL / NA if unavailable
  #   gradient  — score vector at solution (populated by some optimizers); NULL otherwise
  #   vcov      — variance-covariance matrix; NA if Hessian not invertible
  #   counts    — c(fn, gr) evaluation counts from optim()
  #   message   — convergence message string from optim()
  fit_state <- list(
    theta = theta,
    converged = NA,
    objective = NA_real_,
    se = NULL,
    gradient = NULL
  )

  .hzr_run_fit_safely <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("NaNs produced", msg, fixed = TRUE) ||
            grepl("Hessian not invertible; standard errors unavailable", msg, fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  # Distribution dispatch — maximise log-likelihood via the distribution-specific
  # optimiser.  Each branch calls .hzr_optim_<dist>() and writes the standardised
  # fit_state fields.  See the file header for parameterisation conventions and
  # instructions for adding new distributions.
  if (fit && !is.null(theta)) {
    if (dist == "weibull") {
      optim_result <- .hzr_run_fit_safely(.hzr_optim_weibull(
        time = time,
        status = status,
        time_lower = time_lower,
        time_upper = time_upper,
        x = x_fit,
        theta_start = theta,
        control = control
      ))

      fit_state$theta <- optim_result$par
      fit_state$objective <- optim_result$value
      fit_state$converged <- (optim_result$convergence == 0)
      fit_state$se <- if (!anyNA(optim_result$vcov)) sqrt(diag(optim_result$vcov)) else NA
      fit_state$vcov <- optim_result$vcov
      fit_state$counts <- optim_result$counts
      fit_state$message <- optim_result$message
    } else if (dist == "exponential") {
      optim_result <- .hzr_run_fit_safely(.hzr_optim_exponential(
        time = time,
        status = status,
        time_lower = time_lower,
        time_upper = time_upper,
        x = x_fit,
        theta_start = theta,
        control = control
      ))

      fit_state$theta <- optim_result$par
      fit_state$objective <- optim_result$value
      fit_state$converged <- (optim_result$convergence == 0)
      fit_state$se <- if (!anyNA(optim_result$vcov)) sqrt(diag(optim_result$vcov)) else NA
      fit_state$vcov <- optim_result$vcov
      fit_state$counts <- optim_result$counts
      fit_state$message <- optim_result$message
    } else if (dist == "loglogistic") {
      optim_result <- .hzr_run_fit_safely(.hzr_optim_loglogistic(
        time = time,
        status = status,
        time_lower = time_lower,
        time_upper = time_upper,
        x = x_fit,
        theta_start = theta,
        control = control
      ))

      fit_state$theta <- optim_result$par
      fit_state$objective <- optim_result$value
      fit_state$converged <- (optim_result$convergence == 0)
      fit_state$se <- if (!anyNA(optim_result$vcov)) sqrt(diag(optim_result$vcov)) else NA
      fit_state$vcov <- optim_result$vcov
      fit_state$counts <- optim_result$counts
      fit_state$message <- optim_result$message
    } else if (dist == "lognormal") {
      optim_result <- .hzr_run_fit_safely(.hzr_optim_lognormal(
        time = time,
        status = status,
        time_lower = time_lower,
        time_upper = time_upper,
        x = x_fit,
        theta_start = theta,
        control = control
      ))

      fit_state$theta <- optim_result$par
      fit_state$objective <- optim_result$value
      fit_state$converged <- (optim_result$convergence == 0)
      fit_state$se <- if (!anyNA(optim_result$vcov)) sqrt(diag(optim_result$vcov)) else NA
      fit_state$vcov <- optim_result$vcov
      fit_state$counts <- optim_result$counts
      fit_state$message <- optim_result$message
    } else {
      stop("Distribution '", dist, "' not yet supported for fitting.", call. = FALSE)
    }
  }

  # Assemble the hazard S3 object.
  # $call       — captured call for reproducibility / print
  # $spec       — model specification (dist, control)
  # $data       — raw data stored for default predict() / refit
  # $fit        — optimisation results (see fit_state fields above)
  # $legacy_args — pass-through ... args for SAS-migration parity
  # $engine     — implementation tag ("native-r-m2")
  obj <- list(
    call = match.call(),
    spec = list(dist = dist, control = control, time_windows = time_windows),
    data = list(
      time = time,
      time_lower = time_lower,
      time_upper = time_upper,
      status = as.numeric(status),
      x = x
    ),
    fit = fit_state,
    legacy_args = list(...),
    engine = "native-r-m2"
  )

  class(obj) <- "hazard"
  obj
}

#' Predict from a hazard model object
#'
#' Produces prediction outputs from a `hazard` object. Supports multiple prediction
#' types including linear predictor, hazard, survival probability, and cumulative hazard.
#'
#' @param object A `hazard` object.
#' @param newdata Optional matrix or data frame of predictors. For types requiring
#'   time (e.g., "survival", "cumulative_hazard"), newdata should include a `time`
#'   column, or time will be taken from the fitted object's data.
#' @param type Prediction type:
#'   - `"linear_predictor"`: Linear predictor η = x·β
#'   - `"hazard"`: Hazard scale exp(η) (not conditional on time for Weibull)
#'   - `"survival"`: Survival probability S(t|x) = exp(-H(t|x))
#'   - `"cumulative_hazard"`: Cumulative hazard H(t|x) at event times
#' @param ... Unused; included for S3 compatibility.
#'
#' @details
#' For Weibull models with survival or cumulative_hazard predictions:
#' - Cumulative hazard: H(t|x) = (μ·t)^ν · exp(η)
#' - Survival: S(t|x) = exp(-H(t|x))
#'
#' Time values must be positive and finite. If newdata contains a `time` column,
#' it will be used; otherwise, the time vector from the fitted object is used.
#' For models fit with `time_windows`, predictions for `type = "linear_predictor"`
#' or `"hazard"` also require time values (via `newdata$time` or fitted-time fallback)
#' so window-specific coefficients can be selected.
#'
#' @return Numeric vector of predictions.
#' @export
predict.hazard <- function(object, newdata = NULL, type = c("hazard", "linear_predictor", "survival", "cumulative_hazard"), ...) {
  type <- match.arg(type)
  theta <- object$fit$theta
  time_windows <- object$spec$time_windows

  if (is.null(theta)) {
    stop("No coefficients ('theta') are available in 'object'.", call. = FALSE)
  }

  # -----------------------------------------------------------------------
  # Predictions that do NOT need time (linear_predictor, hazard)
  # -----------------------------------------------------------------------
  # These predictions work purely through the linear predictor η = X β.
  # Shape parameters are stripped from theta; only covariate β's are used.
  # hazard returns exp(η), which is the relative-hazard multiplier (not the
  # baseline conditional hazard, which would require time + distribution).
  if (type %in% c("hazard", "linear_predictor")) {
    n_pred <- NULL
    pred_time <- NULL
    if (is.null(newdata)) {
      x <- object$data$x
      pred_time <- object$data$time
    } else {
      newdata <- as.data.frame(newdata)
      n_pred <- nrow(newdata)
      if ("time" %in% names(newdata)) {
        pred_time <- newdata$time
      }
      # Remove time column if present (not needed for hazard/linear_predictor)
      newdata <- newdata[, names(newdata) != "time", drop = FALSE]
      if (ncol(newdata) > 0) {
        x <- as.matrix(newdata)
      } else {
        x <- NULL
      }
    }

    if (!is.null(time_windows)) {
      if (is.null(x)) {
        stop("Time-varying coefficients require predictors for '", type, "' predictions.", call. = FALSE)
      }
      if (is.null(pred_time)) {
        stop("Time-varying coefficients require a 'time' column in newdata for '", type, "' predictions.", call. = FALSE)
      }
      # Apply the same piecewise expansion used at fit time.
      x <- .hzr_expand_time_varying_design(x = x, time = pred_time, time_windows = time_windows)
    }
    
    if (is.null(x)) {
      # No covariates; univariate models with only shape parameters are still valid.
      # Return a zero linear predictor (η = 0 ⇒ no covariate adjustment).
      if (length(theta) == 2 && object$spec$dist == "weibull") {
        # Univariable case: return constant hazard/predictor
        if (is.null(n_pred)) {
          n_pred <- length(object$data$time)
        }
        eta <- rep(0, n_pred)
      } else {
        stop("Predictors are required either in the fitted object or via 'newdata'.", call. = FALSE)
      }
    } else {
      # Extract covariate coefficients from theta
      # For Weibull: theta = [mu, nu, beta1, beta2, ...]
      if (object$spec$dist == "weibull" && length(theta) > 2) {
        # Remove shape parameters, keep only covariate coefficients
        beta <- theta[3:length(theta)]
      } else if (object$spec$dist == "loglogistic" && length(theta) > 2) {
        # Log-logistic: theta = [log(alpha), log(beta), beta1, beta2, ...]
        beta <- theta[3:length(theta)]
      } else if (object$spec$dist == "lognormal" && length(theta) > 2) {
        # Log-normal: theta = [mu, log(sigma), beta1, beta2, ...]
        beta <- theta[3:length(theta)]
      } else if (object$spec$dist == "exponential" && length(theta) > 1) {
        # Exponential: theta = [log(lambda), beta1, beta2, ...]
        beta <- theta[2:length(theta)]
      } else {
        # Assume all of theta are covariate coefficients
        beta <- theta
      }
      
      if (ncol(x) != length(beta)) {
        stop("Number of predictor columns (", ncol(x), 
             ") must match number of covariate coefficients (", length(beta), ").", 
             call. = FALSE)
      }
      eta <- as.numeric(x %*% beta)
    }
    
    if (type == "linear_predictor") {
      return(eta)
    }
    return(exp(eta))
  }

  # -----------------------------------------------------------------------
  # Predictions that DO need time (survival, cumulative_hazard)
  # -----------------------------------------------------------------------
  # Each distribution computes H(t|x) in its own way; all then share the
  # same final return statements:
  #   cumulative_hazard → return(cumhaz)
  #   survival          → return(exp(-cumhaz))
  #
  # EXCEPTION: log-normal uses a direct Φ(−z) formula for survival and
  # returns early inside its branch, bypassing the shared return below.
  # This is because the log-normal is an AFT model where H(t) = −log Φ(−z)
  # is numerically better computed directly from the normal CDF rather than
  # via exp(−H).
  if (type %in% c("survival", "cumulative_hazard")) {
    if (!object$spec$dist %in% c("weibull", "exponential", "loglogistic", "lognormal")) {
      stop("Prediction type '", type, "' is only supported for Weibull, exponential, log-logistic, and log-normal models.", call. = FALSE)
    }

    # Extract time from newdata or use fitted data
    if (!is.null(newdata)) {
      newdata <- as.data.frame(newdata)
      if ("time" %in% names(newdata)) {
        time <- newdata$time
        x <- as.matrix(newdata[, names(newdata) != "time", drop = FALSE])
      } else {
        stop("'newdata' must contain a 'time' column for '", type, "' predictions.", call. = FALSE)
      }
    } else {
      time <- object$data$time
      x <- object$data$x
    }

    if (!is.null(time_windows) && !is.null(x) && ncol(x) > 0) {
      # Keep prediction design matrix consistent with training-time expansion.
      x <- .hzr_expand_time_varying_design(x = x, time = time, time_windows = time_windows)
    }

    if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0)) {
      stop("'time' must be a numeric vector of finite non-negative values.", call. = FALSE)
    }

    # Dispatch by distribution — each branch sets `cumhaz`.
    # Log-normal is the exception: it returns early (see note above).
    if (object$spec$dist == "weibull") {
      # Extract Weibull shape parameters
      if (length(theta) < 2) {
        stop("Weibull model must have at least 2 parameters (mu, nu).", call. = FALSE)
      }

      mu <- theta[1]
      nu <- theta[2]

      if (mu <= 0 || nu <= 0) {
        stop("Weibull shape parameters (mu, nu) must be positive.", call. = FALSE)
      }

      # Compute linear predictor (covariate effect)
      if (!is.null(x) && ncol(x) > 0) {
        if (length(theta) < ncol(x) + 2) {
          stop("Number of parameters insufficient for predictor columns.", call. = FALSE)
        }
        beta <- theta[3:length(theta)]
        eta <- as.numeric(x %*% beta)
      } else {
        eta <- rep(0, length(time))
      }

      # Compute cumulative hazard: H(t|x) = (mu * t)^nu * exp(eta)
      cumhaz <- (mu * time) ^ nu * exp(eta)
    } else if (object$spec$dist == "exponential") {
      # Extract exponential parameter
      if (length(theta) < 1) {
        stop("Exponential model must have at least 1 parameter (log(lambda)).", call. = FALSE)
      }

      log_lambda <- theta[1]
      lambda <- exp(log_lambda)

      # Compute linear predictor (covariate effect)
      if (!is.null(x) && ncol(x) > 0) {
        if (length(theta) < ncol(x) + 1) {
          stop("Number of parameters insufficient for predictor columns.", call. = FALSE)
        }
        beta <- theta[2:length(theta)]
        eta <- as.numeric(x %*% beta)
      } else {
        eta <- rep(0, length(time))
      }

      # Compute cumulative hazard: H(t|x) = lambda * t * exp(eta)
      cumhaz <- lambda * time * exp(eta)
    } else if (object$spec$dist == "loglogistic") {
      # Extract log-logistic parameters
      if (length(theta) < 2) {
        stop("Log-logistic model must have at least 2 parameters (log(alpha), log(beta)).", call. = FALSE)
      }

      log_alpha <- theta[1]
      log_beta <- theta[2]
      alpha <- exp(log_alpha)
      beta <- exp(log_beta)

      # Compute linear predictor (covariate effect)
      if (!is.null(x) && ncol(x) > 0) {
        if (length(theta) < ncol(x) + 2) {
          stop("Number of parameters insufficient for predictor columns.", call. = FALSE)
        }
        beta_coef <- theta[3:length(theta)]
        eta <- as.numeric(x %*% beta_coef)
      } else {
        eta <- rep(0, length(time))
      }

      # Compute cumulative hazard: H(t|x) = log(1 + alpha * t^beta * exp(eta))
      cumhaz <- log(1 + alpha * (time ^ beta) * exp(eta))
    } else if (object$spec$dist == "lognormal") {
      # Extract log-normal parameters (AFT: location + scale)
      if (length(theta) < 2) {
        stop("Log-normal model must have at least 2 parameters (mu, log(sigma)).", call. = FALSE)
      }

      mu <- theta[1]
      log_sigma <- theta[2]
      sigma <- exp(log_sigma)

      # Compute linear predictor (AFT: covariates add to location)
      if (!is.null(x) && ncol(x) > 0) {
        if (length(theta) < ncol(x) + 2) {
          stop("Number of parameters insufficient for predictor columns.", call. = FALSE)
        }
        beta_coef <- theta[3:length(theta)]
        eta <- mu + as.numeric(x %*% beta_coef)
      } else {
        eta <- rep(mu, length(time))
      }

      # Standardized residual: z = (log(t) - eta) / sigma
      z <- (log(time) - eta) / sigma

      # Survival: S(t|x) = Phi(-z)
      surv <- pnorm(-z)

      if (type == "survival") {
        return(surv)
      }

      # Cumulative hazard: H(t|x) = -log(S(t|x))
      return(-pnorm(-z, log.p = TRUE))
    }

    if (type == "cumulative_hazard") {
      return(cumhaz)
    }

    # Survival: S(t|x) = exp(-H(t|x))
    return(exp(-cumhaz))
  }

  stop("Unknown prediction type: '", type, "'.", call. = FALSE)
}


#' @export
print.hazard <- function(x, ...) {
  n <- length(x$data$time)
  p <- if (is.null(x$data$x)) 0L else ncol(x$data$x)
  cat("hazard object\n")
  cat("  observations:", n, "\n")
  cat("  predictors:  ", p, "\n")
  cat("  dist:        ", x$spec$dist, "\n")
  cat("  engine:      ", x$engine, "\n")
  if (!anyNA(x$fit$objective)) {
    cat("  log-lik:     ", format(x$fit$objective, digits = 6), "\n")
    cat("  converged:   ", x$fit$converged, "\n")
  }
  invisible(x)
}

#' Summarize a hazard model
#'
#' Returns a compact summary of a `hazard` object, including model metadata,
#' fit diagnostics, and coefficient-level statistics when available.
#'
#' @param object A `hazard` object.
#' @param ... Unused; for S3 compatibility.
#' @return An object of class `summary.hazard`.
#' @export
summary.hazard <- function(object, ...) {
  n <- length(object$data$time)
  p <- if (is.null(object$data$x)) 0L else ncol(object$data$x)
  theta <- object$fit$theta
  vcov_mat <- object$fit$vcov

  coef_table <- NULL
  if (!is.null(theta)) {
    coef_names <- .hzr_parameter_names(theta = theta, dist = object$spec$dist, p = p)
    std_error <- rep(NA_real_, length(theta))
    z_stat <- rep(NA_real_, length(theta))
    p_value <- rep(NA_real_, length(theta))

    if (!is.null(vcov_mat) && is.matrix(vcov_mat) && !anyNA(vcov_mat)) {
      std_error <- sqrt(diag(vcov_mat))
      valid <- is.finite(std_error) & std_error > 0
      z_stat[valid] <- theta[valid] / std_error[valid]
      p_value[valid] <- 2 * pnorm(-abs(z_stat[valid]))
    }

    coef_table <- data.frame(
      estimate = unname(theta),
      std_error = std_error,
      z_stat = z_stat,
      p_value = p_value,
      row.names = coef_names,
      check.names = FALSE
    )
  }

  out <- list(
    call = object$call,
    n = n,
    p = p,
    dist = object$spec$dist,
    engine = object$engine,
    converged = object$fit$converged,
    log_lik = object$fit$objective,
    counts = object$fit$counts,
    message = object$fit$message,
    coefficients = coef_table,
    has_vcov = !is.null(vcov_mat) && is.matrix(vcov_mat) && !anyNA(vcov_mat)
  )

  class(out) <- "summary.hazard"
  out
}

#' @export
print.summary.hazard <- function(x, ...) {
  cat("hazard model summary\n")
  cat("  observations:", x$n, "\n")
  cat("  predictors:  ", x$p, "\n")
  cat("  dist:        ", x$dist, "\n")
  cat("  engine:      ", x$engine, "\n")

  if (!is.null(x$converged) && !is.na(x$converged)) {
    cat("  converged:   ", x$converged, "\n")
  }
  if (!is.null(x$log_lik) && !is.na(x$log_lik)) {
    cat("  log-lik:     ", format(x$log_lik, digits = 6), "\n")
  }
  if (!is.null(x$counts)) {
    fn_count <- x$counts[["function"]] %||% x$counts[["fn"]] %||% NA_integer_
    gr_count <- x$counts[["gradient"]] %||% x$counts[["gr"]] %||% NA_integer_
    if (!is.na(fn_count) || !is.na(gr_count)) {
      cat("  evaluations: ", "fn=", fn_count, ", gr=", gr_count, "\n", sep = "")
    }
  }
  if (!is.null(x$message) && nzchar(x$message)) {
    cat("  message:     ", x$message, "\n")
  }

  if (!is.null(x$coefficients)) {
    cat("\nCoefficients:\n")
    print(x$coefficients)
  }

  invisible(x)
}

#' Extract coefficients from hazard model
#'
#' @param object A `hazard` object.
#' @param ... Unused; for S3 compatibility.
#' @export
coef.hazard <- function(object, ...) {
  if (is.null(object$fit$theta)) {
    return(NULL)
  }
  object$fit$theta
}

#' Extract variance-covariance matrix from hazard model
#'
#' Returns the estimated variance-covariance matrix of the fitted coefficients.
#'
#' @param object A `hazard` object.
#' @param ... Unused; for S3 compatibility.
#' @export
vcov.hazard <- function(object, ...) {
  if (is.null(object$fit$vcov) || anyNA(object$fit$vcov)) {
    return(NA)
  }
  object$fit$vcov
}

.hzr_parameter_names <- function(theta, dist, p) {
  theta_names <- names(theta)
  if (!is.null(theta_names) && all(nzchar(theta_names))) {
    return(theta_names)
  }

  if (!is.null(p) && p > 0L && length(theta) == p) {
    return(paste0("beta", seq_along(theta)))
  }

  base_names <- switch(
    dist,
    weibull = c("mu", "nu"),
    exponential = c("log_lambda"),
    loglogistic = c("log_alpha", "log_beta"),
    lognormal = c("mu", "log_sigma"),
    character()
  )

  n_theta <- length(theta)
  n_base <- min(length(base_names), n_theta)
  out <- character(n_theta)

  if (n_base > 0) {
    out[seq_len(n_base)] <- base_names[seq_len(n_base)]
  }
  if (n_theta > n_base) {
    out[seq.int(n_base + 1L, n_theta)] <- paste0("beta", seq_len(n_theta - n_base))
  }
  out
}

#' Coerce x to a validated numeric design matrix
#'
#' Accepts data.frame or numeric matrix input; validates dimensions and
#' finiteness.  Called by hazard() to normalise the x argument before any
#' downstream use.
#'
#' @param x A numeric matrix or data frame with numeric columns.
#' @param n Expected row count (length of time/status vectors); checked if non-NULL.
#' @return A numeric matrix with column names preserved (or added if absent).
#' @keywords internal
#'
.hzr_as_design_matrix <- function(x, n = NULL) {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  }


#' Expand predictors for piecewise time-varying coefficients
#'
#' Builds a block design matrix with one predictor block per time window. For each
#' observation, only the block corresponding to that row's time window is active;
#' all other blocks are zero.
#'
#' Example with two predictors and one cut point:
#' - input columns: `x1`, `x2`
#' - windows: `(−Inf, c1]`, `(c1, Inf)`
#' - output columns: `x1_w1`, `x2_w1`, `x1_w2`, `x2_w2`
#'
#' @param x Numeric predictor matrix/data.frame.
#' @param time Numeric time vector used to assign rows to windows.
#' @param time_windows Numeric vector of window cut points.
#' @return Expanded numeric matrix with window-specific column names.
#' @noRd
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("Predictor input must be a numeric matrix or coercible data frame.", call. = FALSE)
  }

  if (any(!is.finite(x))) {
    stop("Predictor matrix contains non-finite values.", call. = FALSE)
  }

  if (!is.null(n) && nrow(x) != n) {
    stop("Predictor rows must match the length of 'time'.", call. = FALSE)
  }

  x
}

.hzr_expand_time_varying_design <- function(x, time, time_windows) {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    stop("Time-varying design requires numeric matrix/data.frame predictors.", call. = FALSE)
  }

  if (!is.numeric(time) || length(time) != nrow(x) || any(!is.finite(time)) || any(time < 0)) {
    stop("'time' must be finite, non-negative, and match predictor rows for time-varying expansion.", call. = FALSE)
  }

  if (!is.numeric(time_windows) || length(time_windows) < 1) {
    stop("'time_windows' must contain at least one cut point.", call. = FALSE)
  }

  time_windows <- sort(unique(time_windows))
  # Bin index k means observation time falls in window k.
  bins <- cut(time, breaks = c(-Inf, time_windows, Inf), right = TRUE, include.lowest = TRUE, labels = FALSE)
  n_bins <- length(time_windows) + 1L

  base_names <- colnames(x)
  if (is.null(base_names)) {
    base_names <- paste0("x", seq_len(ncol(x)))
  }

  out <- matrix(0, nrow = nrow(x), ncol = ncol(x) * n_bins)
  col_ptr <- 1L
  out_names <- character(ncol(out))

  for (k in seq_len(n_bins)) {
    idx <- bins == k
    # Populate only the active window block for rows in this bin.
    block <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    if (any(idx)) {
      block[idx, ] <- x[idx, , drop = FALSE]
    }
    cols <- col_ptr:(col_ptr + ncol(x) - 1L)
    out[, cols] <- block
    out_names[cols] <- paste0(base_names, "_w", k)
    col_ptr <- col_ptr + ncol(x)
  }

  colnames(out) <- out_names
  out
}
