#' Build and optionally fit a hazard model
#'
#' Creates a `hazard` object and optionally fits it via maximum likelihood.
#' This mirrors the argument-oriented workflow of the legacy HAZARD C/SAS
#' implementation: supply starting values in `theta` and the function will
#' optimize to produce fitted estimates.
#'
#' @param time Numeric follow-up time vector.
#' @param status Numeric or logical event indicator vector.
#' @param x Optional design matrix (or data frame coercible to matrix).
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
#' @return An object of class `hazard` with fitted parameters and fit diagnostics.
#' @export
hazard <- function(time,
                   status,
                   x = NULL,
                   theta = NULL,
                   dist = "weibull",
                   fit = FALSE,
                   control = list(),
                   ...) {
  if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0)) {
    stop("'time' must be a numeric vector of finite non-negative values.", call. = FALSE)
  }

  n <- length(time)
  if (length(status) != n) {
    stop("'status' must have the same length as 'time'.", call. = FALSE)
  }

  if (!is.null(x)) {
    x <- .hzr_as_design_matrix(x, n = n)
  }

  if (!is.null(theta)) {
    if (!is.numeric(theta) || any(!is.finite(theta))) {
      stop("'theta' must be a finite numeric vector when provided.", call. = FALSE)
    }

    # If x exists, theta must include coefficients for all variates
    # theta = [shape parms ... | covariate coefficients ...]
    # For now, assume theta length determines whether we expect x
    if (!is.null(x) && length(theta) < ncol(x)) {
      stop("'theta' length must be >= number of columns in 'x'.", call. = FALSE)
    }
  }

  if (!is.character(dist) || length(dist) != 1 || !nzchar(dist)) {
    stop("'dist' must be a non-empty character scalar.", call. = FALSE)
  }

  if (!is.list(control)) {
    stop("'control' must be a list.", call. = FALSE)
  }

  # Initialize fit state (starting values)
  fit_state <- list(
    theta = theta,
    converged = NA,
    objective = NA_real_,
    se = NULL,
    gradient = NULL
  )

  # If fitting requested and theta provided, optimize
  if (fit && !is.null(theta)) {
    if (dist == "weibull") {
      optim_result <- .hzr_optim_weibull(
        time = time,
        status = status,
        x = x,
        theta_start = theta,
        control = control
      )

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

  obj <- list(
    call = match.call(),
    spec = list(dist = dist, control = control),
    data = list(time = time, status = as.numeric(status), x = x),
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
#'
#' @return Numeric vector of predictions.
#' @export
predict.hazard <- function(object, newdata = NULL, type = c("hazard", "linear_predictor", "survival", "cumulative_hazard"), ...) {
  type <- match.arg(type)
  theta <- object$fit$theta

  if (is.null(theta)) {
    stop("No coefficients ('theta') are available in 'object'.", call. = FALSE)
  }

  # For types not requiring time, use simple design matrix approach
  if (type %in% c("hazard", "linear_predictor")) {
    n_pred <- NULL
    if (is.null(newdata)) {
      x <- object$data$x
    } else {
      newdata <- as.data.frame(newdata)
      n_pred <- nrow(newdata)
      # Remove time column if present (not needed for hazard/linear_predictor)
      newdata <- newdata[, names(newdata) != "time", drop = FALSE]
      if (ncol(newdata) > 0) {
        x <- as.matrix(newdata)
      } else {
        x <- NULL
      }
    }
    
    if (is.null(x)) {
      # No covariates; check if theta is just shape parameters
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

  # For types requiring time (survival, cumulative_hazard)
  if (type %in% c("survival", "cumulative_hazard")) {
    if (object$spec$dist != "weibull") {
      stop("Prediction type '", type, "' is only supported for Weibull models.", call. = FALSE)
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

    if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0)) {
      stop("'time' must be a numeric vector of finite non-negative values.", call. = FALSE)
    }

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

.hzr_as_design_matrix <- function(x, n = NULL) {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  }

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
