# math_primitives.R — Numerically stable elementary functions
#
# PURPOSE
# -------
# Floating-point naive implementations of log(1+exp(x)), log(1-exp(-x)), and
# probability clamping produce incorrect results near the boundaries of
# double precision.  These primitives use piece-wise stable formulas that are
# safe across the full range of inputs encountered in survival analysis
# (very small or very large hazard rates, near-zero probabilities).
#
# All three functions are exported so users can call them directly when
# constructing custom likelihoods or diagnostics.

#' Numerically stable log(1 + exp(x))
#'
#' @param x Numeric vector.
#' @return Numeric vector with element-wise `log(1 + exp(x))`.
#' @export
hzr_log1pexp <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

#' Numerically stable log(1 - exp(-x)) for x > 0
#'
#' @param x Numeric vector with positive values.
#' @return Numeric vector with element-wise `log(1 - exp(-x))`.
#' @export
hzr_log1mexp <- function(x) {
  out <- rep(NA_real_, length(x))
  cutoff <- log(2)

  valid <- is.finite(x) & (x > 0)
  small <- valid & (x <= cutoff)
  large <- valid & (x > cutoff)

  out[small] <- log(-expm1(-x[small]))
  out[large] <- log1p(-exp(-x[large]))
  out
}

#' Clamp probabilities away from 0 and 1
#'
#' @param p Numeric vector of probabilities.
#' @param eps Small positive tolerance.
#' @return Numeric vector bounded to `[eps, 1 - eps]`.
#' @export
hzr_clamp_prob <- function(p, eps = 1e-12) {
  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0 || eps >= 0.5) {
    stop("'eps' must be a finite scalar in (0, 0.5).", call. = FALSE)
  }

  pmin(1 - eps, pmax(eps, p))
}
