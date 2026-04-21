#' @keywords internal
NULL

# predict-cl.R -- Delta-method confidence limits for predict.hazard()
#
# PROBLEM
# -------
# predict.hazard() returns a point estimate eta(theta) for each row.  We want
# a standard error and a (1 - alpha) confidence interval around each one,
# matching SAS HAZARD's `hzp_calc_haz_CL.c` / `hzp_calc_srv_CL.c`.
#
# METHOD
# ------
# Delta method: if theta_hat ~ N(theta, V) then any smooth function eta(theta)
# is approximately Gaussian with
#
#   var(eta_hat)  ~  (d eta/d theta)^T  V  (d eta/d theta).
#
# For a batch of n prediction points, we assemble J (n x p) where
# J[i, ] = d eta_i / d theta, then
#
#   se_i  =  sqrt( (J V J^T)_{ii} )  =  sqrt( rowSums( (J V) * J ) ).
#
# SCALE CHOICE
# ------------
# Symmetric CLs on the natural scale can leave a [0, 1] survival or a
# [0, Inf) hazard, so we build them on a transformed scale and back-transform:
#
#   hazard / cumhaz   -> log scale:
#     se(log eta)  =  se(eta) / eta,
#     [eta * exp(-z se_log),  eta * exp(+z se_log)]
#
#   survival          -> log-log scale (log(-log S) = log H):
#     se(log(-log S))  =  se(H) / H,
#     S_lwr  =  exp(-H * exp(+z se_logH)),
#     S_upr  =  exp(-H * exp(-z se_logH))
#
#   linear_predictor  -> natural scale (symmetric):
#     [eta - z se, eta + z se]
#
# ENTRY POINT
# -----------
# `.hzr_predict_with_se()` takes a fitted `hazard` object, a type, and the
# materialised `x` / `time` for prediction.  It dispatches to the appropriate
# distribution-specific jacobian builder, runs the delta-method sandwich, and
# returns a list(fit, se.fit, lower, upper).

# ---------------------------------------------------------------------------
# Jacobian: Weibull
# ---------------------------------------------------------------------------
#
# Parameter layout (natural scale, matching .hzr_logl_weibull):
#   theta = [mu, nu, beta_1, ..., beta_p]
#
# H(t|x) = (mu t)^nu exp(eta),   h_rel = exp(eta),   eta = x beta.
#
# Analytical derivatives used below:
#   dH/dmu       = (nu / mu) H
#   dH/dnu       = log(mu t) H
#   dH/dbeta_j   = x_ij H
#
#   d exp(eta)/dbeta_j = x_ij exp(eta)    (other entries are 0)
#
#   eta itself is linear: d eta / d beta_j = x_ij, rest 0.

#' Jacobian of Weibull predictions with respect to theta
#'
#' @param type Prediction type (see `.hzr_predict_with_se`).
#' @param theta MLE parameter vector `c(mu, nu, beta_1, ...)`.
#' @param time Prediction times (may be NULL for `hazard` / `linear_predictor`).
#' @param x Design matrix (n x p_cov) or NULL.
#' @param p Length of theta.
#' @return Numeric n x p Jacobian.
#' @keywords internal
.hzr_predict_jacobian_weibull <- function(type, theta, time, x, p) {
  n_shape <- 2L
  mu <- theta[1]
  nu <- theta[2]
  beta <- if (p > n_shape) theta[(n_shape + 1L):p] else numeric(0)

  if (!is.null(x) && length(beta) > 0L) {
    eta <- as.numeric(x %*% beta)
  } else if (!is.null(x)) {
    eta <- rep(0, nrow(x))
  } else {
    # No covariates: treat all rows as a single reference
    eta <- if (!is.null(time)) rep(0, length(time)) else 0
  }
  n <- length(eta)
  J <- matrix(0, nrow = n, ncol = p)

  if (type %in% c("cumulative_hazard", "survival")) {
    H <- (mu * time) ^ nu * exp(eta)
    J[, 1L] <- (nu / mu) * H
    J[, 2L] <- log(mu * time) * H
    if (length(beta) > 0L && !is.null(x)) {
      J[, (n_shape + 1L):p] <- x * H
    }
  } else if (type == "hazard") {
    # Current API returns relative hazard = exp(eta); shape params drop out.
    h_rel <- exp(eta)
    if (length(beta) > 0L && !is.null(x)) {
      J[, (n_shape + 1L):p] <- x * h_rel
    }
  } else if (type == "linear_predictor") {
    if (length(beta) > 0L && !is.null(x)) {
      J[, (n_shape + 1L):p] <- x
    }
  }

  J
}

# ---------------------------------------------------------------------------
# Jacobian: multiphase
# ---------------------------------------------------------------------------
#
# H(t|x) = sum_j mu_j(x) Phi_j(t),   mu_j(x) = exp(log_mu_j + x_j beta_j).
#
# For each phase j, its sub-vector in theta contains (in order):
#   log_mu_j, (log_t_half_j, nu_j, m_j) or (log_tau_j, gamma_j, alpha_j, eta_j)
#   or () for a "constant" phase, then beta_j covariates.
#
# Derivatives:
#   dH/d log_mu_j    = mu_j Phi_j
#   dH/d shape_s_j   = mu_j dPhi_j/d shape_s     (from .hzr_phase_derivatives
#                                                 / .hzr_g3_phase_derivatives)
#   dH/d beta_jk     = x_ik mu_j Phi_j

#' Jacobian of multiphase cumulative-hazard predictions
#'
#' Only `cumulative_hazard` and `survival` are delegated here -- `hazard`
#' and `linear_predictor` are rejected upstream for multiphase models.
#'
#' @param theta MLE parameter vector.
#' @param time Prediction times (length n).
#' @param phases Named list of `hzr_phase` objects.
#' @param covariate_counts Named integer vector.
#' @param x_list Named list of per-phase design matrices.
#' @param p Length of theta.
#' @return Numeric n x p Jacobian of H(t|x).
#' @keywords internal
.hzr_predict_jacobian_multiphase <- function(theta, time, phases,
                                              covariate_counts, x_list, p) {
  n <- length(time)
  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  J <- matrix(0, nrow = n, ncol = p)

  pos <- 1L
  for (nm in names(phases)) {
    ph <- phases[[nm]]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], ph)

    # mu_j(x_i)
    if (length(pars$beta) > 0L && !is.null(x_list[[nm]])) {
      eta_j <- pars$log_mu + as.numeric(x_list[[nm]] %*% pars$beta)
    } else {
      eta_j <- rep(pars$log_mu, n)
    }
    mu_j <- exp(eta_j)

    # Phi_j and shape derivatives
    if (ph$type == "constant") {
      Phi_j <- time
      # log_mu only
      J[, pos] <- mu_j * Phi_j
      pos <- pos + 1L
    } else if (ph$type == "g3") {
      tau_j <- exp(pars$log_tau)
      pd <- .hzr_g3_phase_derivatives(time, tau = tau_j,
                                        gamma = pars$gamma,
                                        alpha = pars$alpha,
                                        eta = pars$eta)
      Phi_j <- pd$Phi
      J[, pos] <- mu_j * Phi_j  # log_mu
      J[, pos + 1L] <- mu_j * pd$dPhi_dlog_tau
      J[, pos + 2L] <- mu_j * pd$dPhi_dgamma
      J[, pos + 3L] <- mu_j * pd$dPhi_dalpha
      J[, pos + 4L] <- mu_j * pd$dPhi_deta
      pos <- pos + 5L
    } else {
      t_half_j <- exp(pars$log_t_half)
      pd <- .hzr_phase_derivatives(time, t_half = t_half_j,
                                    nu = pars$nu, m = pars$m,
                                    type = ph$type)
      Phi_j <- pd$Phi
      J[, pos] <- mu_j * Phi_j  # log_mu
      J[, pos + 1L] <- mu_j * pd$dPhi_dlog_thalf
      J[, pos + 2L] <- mu_j * pd$dPhi_dnu
      J[, pos + 3L] <- mu_j * pd$dPhi_dm
      pos <- pos + 4L
    }

    # Covariate betas
    n_beta <- covariate_counts[[nm]]
    if (n_beta > 0L && !is.null(x_list[[nm]])) {
      x_phase <- x_list[[nm]]
      for (k in seq_len(n_beta)) {
        J[, pos] <- x_phase[, k] * mu_j * Phi_j
        pos <- pos + 1L
      }
    } else if (n_beta > 0L) {
      pos <- pos + n_beta
    }
  }

  J
}

# ---------------------------------------------------------------------------
# Numeric jacobian fallback for exp / loglogistic / lognormal
# ---------------------------------------------------------------------------
#
# Rather than hand-code each distribution's closed form, we differentiate the
# prediction function numerically through the existing single-distribution
# cumhaz/survival formulas in predict.hazard().  The prediction function is
# built once per call and returns a length-n vector for any candidate theta.

#' Numeric Jacobian of a single-distribution prediction w.r.t. theta
#'
#' @param predict_fn Function(theta) -> numeric vector of length n.
#' @param theta MLE parameter vector.
#' @return Numeric n x p Jacobian.
#' @keywords internal
.hzr_predict_jacobian_numeric <- function(predict_fn, theta) {
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    return(numDeriv::jacobian(predict_fn, theta))
  }
  eps <- (.Machine$double.eps) ^ (1 / 3)
  p <- length(theta)
  fit <- predict_fn(theta)
  n <- length(fit)
  J <- matrix(0, nrow = n, ncol = p)
  for (j in seq_len(p)) {
    h <- eps * max(abs(theta[j]), 1)
    tp <- theta
    tm <- theta
    tp[j] <- tp[j] + h
    tm[j] <- tm[j] - h
    J[, j] <- (predict_fn(tp) - predict_fn(tm)) / (2 * h)
  }
  J
}

# ---------------------------------------------------------------------------
# Delta-method sandwich + scale transform
# ---------------------------------------------------------------------------

#' Per-row standard errors via the delta method
#'
#' @param J Jacobian (n x p).
#' @param vcov p x p variance-covariance matrix.
#' @return Numeric vector of length n.
#' @keywords internal
.hzr_predict_se_from_jacobian <- function(J, vcov) {
  JV <- J %*% vcov
  se <- sqrt(pmax(rowSums(JV * J), 0))
  se
}

#' Apply a scale transform + back-transform for CLs
#'
#' @param fit Point estimate (numeric vector length n).
#' @param se_nat SE on the natural scale (numeric vector length n).
#' @param level Confidence level.
#' @param scale One of "natural", "log", "loglog_survival".
#' @return data.frame with columns fit, se.fit, lower, upper.
#' @keywords internal
.hzr_predict_cl_from_se <- function(fit, se_nat, level, scale) {
  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha / 2)

  lower <- rep(NA_real_, length(fit))
  upper <- rep(NA_real_, length(fit))

  if (scale == "natural") {
    lower <- fit - z * se_nat
    upper <- fit + z * se_nat
  } else if (scale == "log") {
    # eta > 0 required; any non-positive fit gets NA CLs.
    pos <- is.finite(fit) & fit > 0
    se_log <- ifelse(pos, se_nat / abs(fit), NA_real_)
    lower[pos] <- fit[pos] * exp(-z * se_log[pos])
    upper[pos] <- fit[pos] * exp(z * se_log[pos])
  } else if (scale == "loglog_survival") {
    # fit is S = exp(-H); se_nat is the SE of H (not S).  Build CLs on
    # log(-log S) = log H.
    #   S_lwr = exp(-H * exp(+z * se(log H)))
    #   S_upr = exp(-H * exp(-z * se(log H)))
    # The caller passes `fit = S` and `se_nat = se(H)`; we reconstruct
    # H = -log(S) to compute se(log H) = se(H) / H.
    H <- -log(pmin(pmax(fit, .Machine$double.xmin), 1))
    pos <- is.finite(H) & H > 0
    se_logH <- ifelse(pos, se_nat / H, NA_real_)
    lower[pos] <- exp(-H[pos] * exp(z * se_logH[pos]))
    upper[pos] <- exp(-H[pos] * exp(-z * se_logH[pos]))
  } else {
    stop("Unknown CL scale: '", scale, "'.", call. = FALSE)
  }

  data.frame(fit = fit, se.fit = se_nat, lower = lower, upper = upper)
}

# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

#' Compute point predictions + delta-method CLs
#'
#' Dispatched from predict.hazard() when `se.fit = TRUE`.
#'
#' DELTA-METHOD TARGET
#' -------------------
#' The quantity we actually differentiate depends on the prediction type and
#' the CL scale the type uses:
#'
#'   type = "cumulative_hazard" -> target = H, log-scale CL
#'   type = "survival"          -> target = H, log-log-survival CL
#'                                 (final `fit` is exp(-H))
#'   type = "hazard"            -> target = exp(eta) (single-dist only),
#'                                 log-scale CL
#'   type = "linear_predictor"  -> target = eta, natural-scale CL
#'
#' The caller supplies `diff_fn(theta)` that returns the target vector of
#' length n.  For Weibull and multiphase we build J analytically; for
#' exp / loglogistic / lognormal we fall through to a numeric jacobian of
#' `diff_fn`.
#'
#' @param object Fitted `hazard` object.
#' @param type One of "hazard", "linear_predictor", "survival",
#'   "cumulative_hazard".
#' @param time Numeric prediction times or NULL (for type = "hazard" /
#'   "linear_predictor").
#' @param x Design matrix (single-distribution) or NULL.  Ignored for
#'   multiphase; use `x_list` instead.
#' @param x_list Named list of per-phase design matrices (multiphase only).
#' @param cov_counts Named integer covariate counts (multiphase only).
#' @param phases Named list of `hzr_phase` objects (multiphase only).
#' @param level Confidence level (default 0.95).
#' @param diff_fn Required function(theta) -> numeric vector of length n
#'   returning the delta-method target (H, exp(eta), or eta depending on
#'   `type`).  Used for the point estimate AND for the numeric jacobian
#'   fallback in exp / loglogistic / lognormal.
#' @return data.frame with columns `fit`, `se.fit`, `lower`, `upper`.
#' @keywords internal
.hzr_predict_with_se <- function(object, type, time = NULL,
                                   x = NULL, x_list = NULL,
                                   cov_counts = NULL, phases = NULL,
                                   level = 0.95, diff_fn) {

  theta <- object$fit$theta
  p <- length(theta)
  dist <- object$spec$dist
  target <- diff_fn(theta)

  vcov_mat <- object$fit$vcov
  vcov_ok <- !is.null(vcov_mat) && is.matrix(vcov_mat) &&
               nrow(vcov_mat) == p && ncol(vcov_mat) == p
  if (!vcov_ok) {
    warning("Variance-covariance matrix is unavailable; ",
            "standard errors and CLs will be NA.", call. = FALSE)
    n <- length(target)
    fit <- if (type == "survival") exp(-target) else target
    na_vec <- rep(NA_real_, n)
    return(data.frame(fit = fit, se.fit = na_vec,
                      lower = na_vec, upper = na_vec))
  }

  # Fixed parameters (from `fixed = "shapes"` or CoE) leave NA rows/cols in
  # the expanded vcov.  Treat them as known-with-zero-variance: restrict the
  # delta-method sandwich to the free submatrix of vcov and the
  # corresponding columns of J.
  diag_vcov <- diag(vcov_mat)
  free_idx <- which(is.finite(diag_vcov))
  if (length(free_idx) < p) {
    # At least one parameter is fixed -- check that the off-diagonals among
    # free params are finite; if they aren't, something is wrong and we
    # fall back to NA CLs.
    free_submat <- vcov_mat[free_idx, free_idx, drop = FALSE]
    if (anyNA(free_submat)) {
      warning("Variance-covariance matrix has NA entries outside the ",
              "fixed-parameter rows/cols; standard errors and CLs will be NA.",
              call. = FALSE)
      n <- length(target)
      fit <- if (type == "survival") exp(-target) else target
      na_vec <- rep(NA_real_, n)
      return(data.frame(fit = fit, se.fit = na_vec,
                        lower = na_vec, upper = na_vec))
    }
    vcov_use <- free_submat
  } else if (anyNA(vcov_mat)) {
    warning("Variance-covariance matrix has NA entries; ",
            "standard errors and CLs will be NA.", call. = FALSE)
    n <- length(target)
    fit <- if (type == "survival") exp(-target) else target
    na_vec <- rep(NA_real_, n)
    return(data.frame(fit = fit, se.fit = na_vec,
                      lower = na_vec, upper = na_vec))
  } else {
    vcov_use <- vcov_mat
  }

  # --- Build Jacobian of the delta-method target ---------------------------
  if (dist == "weibull") {
    J <- .hzr_predict_jacobian_weibull(type, theta, time, x, p)
  } else if (dist == "multiphase") {
    # hazard / linear_predictor are rejected upstream for multiphase.
    J <- .hzr_predict_jacobian_multiphase(theta, time, phases,
                                            cov_counts, x_list, p)
  } else {
    J <- .hzr_predict_jacobian_numeric(diff_fn, theta)
  }

  # Restrict J to the free columns so the sandwich dimensions match.
  J <- J[, free_idx, drop = FALSE]

  # --- Sandwich -------------------------------------------------------------
  se_target <- .hzr_predict_se_from_jacobian(J, vcov_use)

  # --- Scale transform -----------------------------------------------------
  if (type == "survival") {
    fit <- exp(-target)
    .hzr_predict_cl_from_se(fit, se_target, level, "loglog_survival")
  } else if (type == "cumulative_hazard" || type == "hazard") {
    .hzr_predict_cl_from_se(target, se_target, level, "log")
  } else if (type == "linear_predictor") {
    .hzr_predict_cl_from_se(target, se_target, level, "natural")
  } else {
    stop("Unknown prediction type '", type, "' for se.fit.", call. = FALSE)
  }
}
