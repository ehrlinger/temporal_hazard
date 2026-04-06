#' Build a hazard model object
#'
#' Creates an initial `hazard` object designed to mirror the argument-oriented
#' workflow used by the legacy HAZARD C/SAS implementation. This first version
#' stores validated inputs and a minimal fit state that downstream methods can
#' consume as additional ported logic is implemented.
#'
#' @param time Numeric follow-up time vector.
#' @param status Numeric or logical event indicator vector.
#' @param x Optional design matrix (or data frame coercible to matrix).
#' @param theta Optional numeric coefficient vector.
#' @param dist Character baseline distribution label.
#' @param control Named list of control options.
#' @param ... Additional named arguments retained for parity with legacy calling
#'   conventions.
#'
#' @return An object of class `hazard`.
#' @export
hazard <- function(time,
                   status,
                   x = NULL,
                   theta = NULL,
                   dist = "weibull",
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

    if (!is.null(x) && length(theta) != ncol(x)) {
      stop("'theta' length must match the number of columns in 'x'.", call. = FALSE)
    }
  }

  if (!is.character(dist) || length(dist) != 1 || !nzchar(dist)) {
    stop("'dist' must be a non-empty character scalar.", call. = FALSE)
  }

  if (!is.list(control)) {
    stop("'control' must be a list.", call. = FALSE)
  }

  obj <- list(
    call = match.call(),
    spec = list(dist = dist, control = control),
    data = list(time = time, status = as.numeric(status), x = x),
    fit = list(theta = theta, converged = NA, objective = NA_real_),
    legacy_args = list(...),
    engine = "native-r-stub"
  )

  class(obj) <- "hazard"
  obj
}

#' Predict from a hazard model object
#'
#' Produces basic prediction outputs from a `hazard` object. This is an initial
#' compatibility method that will evolve toward full `hazpred` parity.
#'
#' @param object A `hazard` object.
#' @param newdata Optional matrix or data frame of predictors.
#' @param type Prediction type: `"hazard"` or `"linear_predictor"`.
#' @param ... Unused; included for S3 compatibility.
#'
#' @return Numeric vector of predictions.
#' @export
predict.hazard <- function(object, newdata = NULL, type = c("hazard", "linear_predictor"), ...) {
  type <- match.arg(type)
  theta <- object$fit$theta

  if (is.null(theta)) {
    stop("No coefficients ('theta') are available in 'object'.", call. = FALSE)
  }

  x <- if (is.null(newdata)) object$data$x else .hzr_as_design_matrix(newdata)
  if (is.null(x)) {
    stop("Predictors are required either in the fitted object or via 'newdata'.", call. = FALSE)
  }

  if (ncol(x) != length(theta)) {
    stop("Predictor columns must match coefficient length.", call. = FALSE)
  }

  eta <- as.numeric(x %*% theta)
  if (type == "linear_predictor") {
    return(eta)
  }

  exp(eta)
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
  invisible(x)
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
