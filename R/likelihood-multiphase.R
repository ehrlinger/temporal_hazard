# likelihood-multiphase.R — Multiphase additive hazard likelihood and optimizer
#
# REFERENCES
# ----------
# Blackstone EH, Naftel DC, Turner ME Jr. The decomposition of time-varying
# hazard into phases, each incorporating a separate stream of concomitant
# information. J Am Stat Assoc. 1986;81(395):615-624.
#
# Rajeswaran J, Blackstone EH, Ehrlinger J, Li L, Ishwaran H, Parides MK.
# Probability of atrial fibrillation after ablation: Using a parametric
# nonlinear temporal decomposition mixed effects model. Stat Methods Med Res.
# 2018;27(1):126-141.
#
# MODEL
# -----
# Additive cumulative hazard across J phases:
#
#   H(t | x) = sum_{j=1}^{J}  mu_j(x) * Phi_j(t; t_half_j, nu_j, m_j)
#
# where:
#   mu_j(x) = exp(alpha_j + x_j * beta_j)  — phase-specific log-linear scale
#   Phi_j(t) depends on phase type:
#     "cdf"      → G(t)             bounded [0, 1]   (early risk)
#     "hazard"   → -log(1 - G(t))   monotone incr.   (late risk)
#     "constant" → t                flat rate         (background)
#
# Total hazard:
#   h(t | x) = sum_j  mu_j(x) * phi_j(t)
#
# where phi_j = dPhi_j/dt.
#
# Survival:
#   S(t | x) = exp(-H(t | x))
#
# THETA LAYOUT (internal / estimation scale)
# ------------------------------------------
# For each phase j, the sub-vector is:
#   [log_mu_j, log_t_half_j, nu_j, m_j, beta_j_1, ..., beta_j_pj]
#   (constant phases omit log_t_half, nu, m)
#
# The full theta is the concatenation of all phase sub-vectors.
#
# FUNCTIONS
# ---------
#   .hzr_logl_multiphase()     — log-likelihood
#   .hzr_optim_multiphase()    — optimizer using .hzr_optim_generic()
#   .hzr_multiphase_cumhaz()   — compute total H(t|x) from theta + phases
#   .hzr_multiphase_hazard()   — compute total h(t|x) from theta + phases
#   .hzr_split_theta()         — split theta into per-phase sub-vectors


# ============================================================================
# Theta vector manipulation
# ============================================================================

#' Split the full theta vector into per-phase sub-vectors
#'
#' @param theta Numeric vector — full parameter vector (internal scale).
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector — number of covariates per phase.
#' @return Named list of numeric vectors, one per phase.
#' @keywords internal
.hzr_split_theta <- function(theta, phases, covariate_counts) {
  result <- vector("list", length(phases))
  names(result) <- names(phases)
  pos <- 1L

  for (nm in names(phases)) {
    np <- .hzr_phase_n_params(phases[[nm]], n_covariates = covariate_counts[[nm]])
    result[[nm]] <- theta[pos:(pos + np - 1L)]
    pos <- pos + np
  }

  result
}


#' Extract parameters from a phase sub-vector (internal scale)
#'
#' @param theta_j Numeric sub-vector for one phase.
#' @param phase An `hzr_phase` object.
#' @return Named list with `log_mu`, `log_t_half`, `nu`, `m`, `beta`.
#' @keywords internal
.hzr_unpack_phase_theta <- function(theta_j, phase) {
  pos <- 1L
  log_mu <- theta_j[pos]
  pos <- pos + 1L

  if (phase$type == "g3") {
    log_tau <- theta_j[pos]
    pos <- pos + 1L
    gamma_  <- theta_j[pos]
    pos <- pos + 1L
    alpha_  <- theta_j[pos]
    pos <- pos + 1L
    eta_    <- theta_j[pos]
    pos <- pos + 1L
    log_t_half <- NA_real_
    nu <- NA_real_
    m <- NA_real_
  } else if (phase$type != "constant") {
    log_t_half <- theta_j[pos]
    pos <- pos + 1L
    nu <- theta_j[pos]
    pos <- pos + 1L
    m <- theta_j[pos]
    pos <- pos + 1L
    log_tau <- NA_real_
    gamma_ <- NA_real_
    alpha_ <- NA_real_
    eta_ <- NA_real_
  } else {
    log_t_half <- NA_real_
    nu <- NA_real_
    m <- NA_real_
    log_tau <- NA_real_
    gamma_ <- NA_real_
    alpha_ <- NA_real_
    eta_ <- NA_real_
  }

  beta <- if (pos <= length(theta_j)) theta_j[pos:length(theta_j)] else numeric(0)

  list(log_mu = log_mu, log_t_half = log_t_half, nu = nu, m = m,
       log_tau = log_tau, gamma = gamma_, alpha = alpha_, eta = eta_,
       beta = beta)
}


# ============================================================================
# Cumulative hazard and hazard computations
# ============================================================================

#' Compute total cumulative hazard H(t|x) for multiphase model
#'
#' @param time Numeric vector of times (n).
#' @param theta Full parameter vector (internal scale).
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector of covariate counts per phase.
#' @param x_list Named list of design matrices, one per phase. Each element is
#'   an n x p_j matrix or NULL.
#' @param per_phase Logical; if TRUE return a list with total and per-phase
#'   contributions.
#' @return If `per_phase = FALSE`: numeric vector of length n (total H(t|x)).
#'   If `per_phase = TRUE`: named list with `$total` and one element per phase.
#' @keywords internal
.hzr_multiphase_cumhaz <- function(time, theta, phases, covariate_counts,
                                    x_list, per_phase = FALSE) {
  n <- length(time)
  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  total <- rep(0, n)
  phase_contributions <- if (per_phase) vector("list", length(phases)) else NULL

  for (i in seq_along(phases)) {
    nm <- names(phases)[i]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])

    # mu_j(x) = exp(log_mu + x_j * beta_j)
    if (length(pars$beta) > 0 && !is.null(x_list[[nm]])) {
      eta_j <- pars$log_mu + as.numeric(x_list[[nm]] %*% pars$beta)
    } else {
      eta_j <- rep(pars$log_mu, n)
    }
    mu_j <- exp(eta_j)

    # Phi_j(t)
    if (phases[[nm]]$type == "constant") {
      phi_j <- time
    } else if (phases[[nm]]$type == "g3") {
      tau_j <- exp(pars$log_tau)
      d3 <- hzr_decompos_g3(time, tau = tau_j, gamma = pars$gamma,
                              alpha = pars$alpha, eta = pars$eta)
      phi_j <- d3$G3
    } else {
      t_half_j <- exp(pars$log_t_half)
      # Guard against invalid decomposition parameters
      if (pars$m < 0 && pars$nu < 0) return(rep(Inf, n))
      phi_j <- hzr_phase_cumhaz(time, t_half = t_half_j, nu = pars$nu,
                                  m = pars$m, type = phases[[nm]]$type)
    }

    contrib <- mu_j * phi_j
    total <- total + contrib

    if (per_phase) phase_contributions[[i]] <- contrib
  }

  if (per_phase) {
    names(phase_contributions) <- names(phases)
    phase_contributions$total <- total
    return(phase_contributions)
  }

  total
}


#' Compute total instantaneous hazard h(t|x) for multiphase model
#'
#' @inheritParams .hzr_multiphase_cumhaz
#' @return Numeric vector of length n.
#' @keywords internal
.hzr_multiphase_hazard <- function(time, theta, phases, covariate_counts,
                                    x_list) {
  n <- length(time)
  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  total <- rep(0, n)

  for (i in seq_along(phases)) {
    nm <- names(phases)[i]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])

    if (length(pars$beta) > 0 && !is.null(x_list[[nm]])) {
      eta_j <- pars$log_mu + as.numeric(x_list[[nm]] %*% pars$beta)
    } else {
      eta_j <- rep(pars$log_mu, n)
    }
    mu_j <- exp(eta_j)

    if (phases[[nm]]$type == "constant") {
      dphi_j <- rep(1, n)
    } else if (phases[[nm]]$type == "g3") {
      tau_j <- exp(pars$log_tau)
      d3 <- hzr_decompos_g3(time, tau = tau_j, gamma = pars$gamma,
                              alpha = pars$alpha, eta = pars$eta)
      dphi_j <- d3$g3
    } else {
      t_half_j <- exp(pars$log_t_half)
      if (pars$m < 0 && pars$nu < 0) return(rep(Inf, n))
      dphi_j <- hzr_phase_hazard(time, t_half = t_half_j, nu = pars$nu,
                                   m = pars$m, type = phases[[nm]]$type)
    }

    total <- total + mu_j * dphi_j
  }

  total
}


# ============================================================================
# Conservation of Events (CoE) — Turner's theorem
# ============================================================================
#
# For the additive hazard model H(t|x) = Sigma_j mu_j(x) * Phi_j(t),
# the MLE constraint (score equation for mu) implies:
#
#   Sigma_i E(i) = Sigma_i H(t_i | x_i)
#
# i.e., total observed events equals total predicted cumulative hazard.
#
# This allows one phase's log_mu to be solved analytically at each
# optimizer iteration, reducing the optimization dimension by 1.
#
# Reference: Turner ME Jr., as implemented in C HAZARD setcoe.c / consrv.c.

#' Find the position of each phase's log_mu in the full theta vector
#'
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector of per-phase covariate counts.
#' @return Named integer vector: position of log_mu for each phase.
#' @keywords internal
.hzr_log_mu_positions <- function(phases, covariate_counts) {
  positions <- setNames(integer(length(phases)), names(phases))
  pos <- 1L
  for (nm in names(phases)) {
    positions[[nm]] <- pos
    np <- .hzr_phase_n_params(phases[[nm]],
                               n_covariates = covariate_counts[[nm]])
    pos <- pos + np
  }
  positions
}


#' Select the phase whose log_mu will be solved by conservation
#'
#' Chooses the phase contributing the largest share of total cumulative
#' hazard, matching the C HAZARD SETCOE strategy.
#'
#' @param theta Full parameter vector (internal scale).
#' @param time Numeric vector of follow-up times.
#' @param status Numeric event indicator.
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector.
#' @param x_list Named list of per-phase design matrices.
#' @return Character: name of the phase to fix.
#' @keywords internal
.hzr_select_fixmu_phase <- function(theta, time, status,
                                     phases, covariate_counts, x_list) {
  decomp <- .hzr_multiphase_cumhaz(time, theta, phases,
                                     covariate_counts, x_list,
                                     per_phase = TRUE)
  phase_sums <- vapply(names(phases), function(nm) {
    sum(decomp[[nm]])
  }, numeric(1))

  names(which.max(phase_sums))
}


#' Apply the Conservation of Events adjustment to one phase's log_mu
#'
#' Given the current theta vector, analytically solve the fixmu phase's
#' log_mu so that total predicted events = total observed events.
#'
#' This is called BEFORE each likelihood evaluation inside the optimizer,
#' matching the C HAZARD `CONSRV` entry point.
#'
#' @param theta Full parameter vector (internal scale).
#' @param fixmu_phase Character: name of the phase whose log_mu is solved.
#' @param fixmu_pos Integer: position of that log_mu in theta.
#' @param time Numeric vector of follow-up times.
#' @param status Numeric event indicator.
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector.
#' @param x_list Named list of per-phase design matrices.
#' @param total_events Numeric: sum of observed events (precomputed).
#' @return Updated theta vector with fixmu phase's log_mu adjusted.
#' @keywords internal
.hzr_conserve_events <- function(theta, fixmu_phase, fixmu_pos,
                                  time, status,
                                  phases, covariate_counts, x_list,
                                  total_events) {
  # Compute per-phase cumulative hazard contributions
  decomp <- .hzr_multiphase_cumhaz(time, theta, phases,
                                     covariate_counts, x_list,
                                     per_phase = TRUE)

  # Total predicted events (sum of cumhaz across all observations)
  sumcz <- sum(decomp$total)

  # Contribution from the fixmu phase alone
  sumcj <- sum(decomp[[fixmu_phase]])

  if (sumcj <= 0 || !is.finite(sumcj)) return(theta)

  # Discrepancy: how many events are unaccounted for
  devent <- total_events - sumcz

  # Events the fixmu phase should absorb
  jevent <- sumcj + devent

  if (jevent <= 0 || !is.finite(jevent)) return(theta)

  # Multiplicative adjustment in log scale
  lfactor <- log(jevent / sumcj)

  if (!is.finite(lfactor)) return(theta)

  theta[fixmu_pos] <- theta[fixmu_pos] + lfactor
  theta
}


# ============================================================================
# Log-likelihood
# ============================================================================

#' Log-likelihood for multiphase additive hazard model
#'
#' @param theta Full parameter vector (internal scale).
#' @param time Numeric vector of follow-up times (n).
#' @param status Numeric event indicator: 1 = event, 0 = right-censored,
#'   -1 = left-censored, 2 = interval-censored.
#' @param time_lower Optional lower bounds for interval censoring.
#' @param time_upper Optional upper bounds for left/interval censoring.
#' @param x Design matrix (unused directly; kept for interface compatibility).
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector of per-phase covariate counts.
#' @param x_list Named list of per-phase design matrices.
#' @param return_gradient Logical (ignored; gradient via separate function).
#' @param return_hessian Logical (ignored).
#' @param ... Ignored.
#' @return Scalar log-likelihood. Returns -Inf for infeasible parameters.
#' @keywords internal
.hzr_logl_multiphase <- function(theta, time, status,
                                  time_lower = NULL, time_upper = NULL,
                                  x = NULL,
                                  phases, covariate_counts, x_list,
                                  return_gradient = FALSE,
                                  return_hessian = FALSE, ...) {

  # Feasibility: check parameter constraints

  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  for (nm in names(phases)) {
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])
    if (phases[[nm]]$type %in% c("cdf", "hazard")) {
      if (pars$m < 0 && pars$nu < 0) return(-Inf)
    }
    if (phases[[nm]]$type == "g3") {
      # G3 requires gamma > 0, alpha >= 0, eta > 0
      if (pars$gamma <= 0 || pars$eta <= 0 || pars$alpha < 0) return(-Inf)
    }
  }

  # Validate event times
  if (any(status == 1 & time <= 0)) return(-Inf)

  # Normalize censoring bounds
  lower <- if (is.null(time_lower)) time else time_lower
  upper <- if (is.null(time_upper)) time else time_upper

  # Cumulative hazard at event/censoring times
  cumhaz <- .hzr_multiphase_cumhaz(time, theta, phases, covariate_counts, x_list)
  if (any(!is.finite(cumhaz))) return(-Inf)

  # Instantaneous hazard at event times (only needed for exact events)
  idx_event <- status == 1

  logl <- 0

  # Exact events: log h(t) - H(t)
  if (any(idx_event)) {
    haz <- .hzr_multiphase_hazard(time, theta, phases, covariate_counts, x_list)
    if (any(!is.finite(haz[idx_event])) || any(haz[idx_event] <= 0)) return(-Inf)
    logl <- logl + sum(log(haz[idx_event])) - sum(cumhaz[idx_event])
  }

  # Right-censored: -H(t)
  idx_right <- status == 0
  if (any(idx_right)) {
    logl <- logl - sum(cumhaz[idx_right])
  }

  # Left-censored: log(1 - exp(-H(u)))
  idx_left <- status == -1
  if (any(idx_left)) {
    cumhaz_upper <- .hzr_multiphase_cumhaz(upper, theta, phases,
                                             covariate_counts, x_list)
    logl <- logl + sum(hzr_log1mexp(cumhaz_upper[idx_left]))
  }

  # Interval-censored: log(S(l) - S(u)) = -H(l) + log(1 - exp(-(H(u) - H(l))))
  idx_interval <- status == 2
  if (any(idx_interval)) {
    cumhaz_lower <- .hzr_multiphase_cumhaz(lower, theta, phases,
                                             covariate_counts, x_list)
    cumhaz_upper_iv <- .hzr_multiphase_cumhaz(upper, theta, phases,
                                                covariate_counts, x_list)
    delta_h <- cumhaz_upper_iv[idx_interval] - cumhaz_lower[idx_interval]
    logl <- logl + sum(-cumhaz_lower[idx_interval] + hzr_log1mexp(delta_h))
  }

  if (!is.finite(logl)) return(-Inf)

  logl
}


# ============================================================================
# Analytic gradient
# ============================================================================

#' Gradient of the multiphase log-likelihood
#'
#' Computes the score vector \eqn{d\ell / d\theta} using analytic chain-rule
#' formulas for `log_mu` and `beta` parameters, and central-difference
#' derivatives (via `.hzr_phase_derivatives()`) for shape parameters
#' (`log_t_half`, `nu`, `m`).
#'
#' @inheritParams .hzr_logl_multiphase
#' @return Numeric vector of length `length(theta)` — the gradient.
#'   Returns a zero vector if any component is non-finite (guards optimizer).
#' @keywords internal
.hzr_gradient_multiphase <- function(theta, time, status,
                                      time_lower = NULL, time_upper = NULL,
                                      x = NULL,
                                      phases, covariate_counts, x_list,
                                      ...) {
  n <- length(time)
  p <- length(theta)
  grad <- numeric(p)

  # Feasibility check
  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  for (nm in names(phases)) {
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])
    if (phases[[nm]]$type %in% c("cdf", "hazard")) {
      if (pars$m < 0 && pars$nu < 0) return(grad)
    }
    if (phases[[nm]]$type == "g3") {
      if (pars$gamma <= 0 || pars$eta <= 0 || pars$alpha < 0) return(grad)
    }
  }

  # Normalize censoring bounds
  lower <- if (is.null(time_lower)) time else time_lower
  upper <- if (is.null(time_upper)) time else time_upper

  # Observation type masks
  idx_event    <- which(status == 1)
  idx_right    <- which(status == 0)
  idx_left     <- which(status == -1)

  idx_interval <- which(status == 2)

  has_left     <- length(idx_left) > 0
  has_interval <- length(idx_interval) > 0

  # ── Per-phase quantities ──────────────────────────────────────────────────
  # We need, for each phase j and each observation i:
  #   mu_j(x_i), Phi_j(t_i), phi_j(t_i), and shape derivatives of Phi and phi
  #
  # Then the per-observation gradient contributions are assembled from:
  #   H(t) = sum_j mu_j * Phi_j    (cumulative hazard)
  #   h(t) = sum_j mu_j * phi_j    (instantaneous hazard)

  # Pre-compute total H(t) and h(t) at observation times (and censoring bounds)
  H_t <- rep(0, n)      # total cumulative hazard at time[i]
  h_t <- rep(0, n)      # total instantaneous hazard at time[i]

  # Storage for per-phase intermediate quantities
  n_phases <- length(phases)
  phase_mu    <- vector("list", n_phases)  # mu_j(x_i) for each i
  phase_Phi   <- vector("list", n_phases)  # Phi_j(t_i) for each i
  phase_phi   <- vector("list", n_phases)  # phi_j(t_i) for each i
  phase_deriv <- vector("list", n_phases)  # full derivative list from .hzr_phase_derivatives

  for (j in seq_along(phases)) {
    nm <- names(phases)[j]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])

    # mu_j(x_i) = exp(log_mu + x * beta)
    if (length(pars$beta) > 0 && !is.null(x_list[[nm]])) {
      eta_j <- pars$log_mu + as.numeric(x_list[[nm]] %*% pars$beta)
    } else {
      eta_j <- rep(pars$log_mu, n)
    }
    mu_j <- exp(eta_j)
    phase_mu[[j]] <- mu_j

    # Phi_j, phi_j, and shape derivatives
    if (phases[[nm]]$type == "constant") {
      phase_Phi[[j]] <- time
      phase_phi[[j]] <- rep(1, n)
      phase_deriv[[j]] <- NULL  # no shape params for constant
    } else if (phases[[nm]]$type == "g3") {
      tau_j <- exp(pars$log_tau)
      d3 <- hzr_decompos_g3(time, tau = tau_j, gamma = pars$gamma,
                              alpha = pars$alpha, eta = pars$eta)
      phase_Phi[[j]] <- d3$G3
      phase_phi[[j]] <- d3$g3
      # G3 shape derivatives via finite differences (for now)
      phase_deriv[[j]] <- .hzr_g3_phase_derivatives(
        time, tau = tau_j, gamma = pars$gamma,
        alpha = pars$alpha, eta = pars$eta)
    } else {
      t_half_j <- exp(pars$log_t_half)
      pd <- .hzr_phase_derivatives(time, t_half = t_half_j, nu = pars$nu,
                                    m = pars$m, type = phases[[nm]]$type)
      phase_Phi[[j]]   <- pd$Phi
      phase_phi[[j]]   <- pd$phi
      phase_deriv[[j]] <- pd
    }

    H_t <- H_t + mu_j * phase_Phi[[j]]
    h_t <- h_t + mu_j * phase_phi[[j]]
  }

  # Guard: if total hazard or cumhaz is non-finite, return zero gradient
  if (any(!is.finite(H_t)) || any(!is.finite(h_t))) return(grad)

  # ── Per-observation weight vectors ────────────────────────────────────────
  # dLogl/dH(t_i) depends on observation type:
  #   event:         d/dH [ log h(t) - H(t) ] = -1        (H part)
  #   right-cens:    d/dH [ -H(t) ]            = -1
  #   left-cens:     d/dH [ log(1-exp(-H(u))) ] = exp(-H) / (1 - exp(-H))
  #   interval-cens: more complex (handled separately)
  #
  # For events, there's also the d/d[theta] of log h(t):
  #   d(log h) / d(theta_j) = (1/h) * dh/d(theta_j)

  # Weight for the cumulative hazard part: w_H_i = dLogl / dH(t_i)
  w_H <- numeric(n)
  w_H[idx_event] <- -1
  w_H[idx_right] <- -1

  if (has_left) {
    H_upper <- .hzr_multiphase_cumhaz(upper, theta, phases, covariate_counts,
                                       x_list)
    # d/dH log(1 - exp(-H)) = exp(-H) / (1 - exp(-H))
    exp_neg_H_u <- exp(-H_upper[idx_left])
    one_minus   <- pmax(1 - exp_neg_H_u, .Machine$double.xmin)
    w_H[idx_left] <- exp_neg_H_u / one_minus
    # Left-censored uses H(upper), not H(time), so we handle below
  }

  # Weight for the instantaneous hazard part (events only): 1/h(t_i)
  # This multiplies dh/d(theta_j)
  inv_h <- numeric(n)
  if (length(idx_event) > 0) {
    h_event <- h_t[idx_event]
    h_event <- pmax(h_event, .Machine$double.xmin)
    inv_h[idx_event] <- 1 / h_event
  }

  # ── Assemble gradient per phase ───────────────────────────────────────────
  pos <- 1L  # position in the full theta vector

  for (j in seq_along(phases)) {
    nm <- names(phases)[j]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])
    mu_j  <- phase_mu[[j]]
    Phi_j <- phase_Phi[[j]]
    phi_j <- phase_phi[[j]]

    # ── d(logl) / d(log_mu_j) ──────────────────────────────────────────
    # From H part: w_H_i * dH/d(log_mu) = w_H_i * mu_j * Phi_j
    # From h part (events): inv_h_i * dh/d(log_mu) = inv_h_i * mu_j * phi_j
    dlogl_dlog_mu <- sum(w_H * mu_j * Phi_j) +
                     sum(inv_h[idx_event] * mu_j[idx_event] * phi_j[idx_event])

    # Left-censored: uses H(upper), need separate Phi_j evaluation there
    if (has_left) {
      if (phases[[nm]]$type == "constant") {
        Phi_j_upper <- upper[idx_left]
      } else if (phases[[nm]]$type == "g3") {
        tau_j <- exp(pars$log_tau)
        d3_upper <- hzr_decompos_g3(upper[idx_left], tau = tau_j,
                                      gamma = pars$gamma, alpha = pars$alpha,
                                      eta = pars$eta)
        Phi_j_upper <- d3_upper$G3
      } else {
        t_half_j <- exp(pars$log_t_half)
        Phi_j_upper <- hzr_phase_cumhaz(upper[idx_left], t_half = t_half_j,
                                          nu = pars$nu, m = pars$m,
                                          type = phases[[nm]]$type)
      }
      dlogl_dlog_mu <- dlogl_dlog_mu +
        sum(w_H[idx_left] * mu_j[idx_left] * Phi_j_upper)
      # Subtract the w_H * mu * Phi at time[i] that we already added above
      # (for left-censored, the H contribution is at upper, not time)
      dlogl_dlog_mu <- dlogl_dlog_mu -
        sum(w_H[idx_left] * mu_j[idx_left] * Phi_j[idx_left])
    }

    # Interval-censored: use numerical fallback for these observations
    # (handled at end of loop via finite-difference correction)

    grad[pos] <- dlogl_dlog_mu
    pos <- pos + 1L

    # ── Shape parameters (non-constant phases only) ─────────────────────
    if (phases[[nm]]$type %in% c("cdf", "hazard")) {
      pd <- phase_deriv[[j]]

      # Compute shape derivatives at upper[idx_left] for left-censoring
      pd_upper <- NULL
      if (has_left) {
        t_half_j <- exp(pars$log_t_half)
        pd_upper <- .hzr_phase_derivatives(upper[idx_left],
          t_half = t_half_j, nu = pars$nu, m = pars$m,
          type = phases[[nm]]$type)
      }

      # For each shape param s in {log_t_half, nu, m}:
      shape_derivs <- list(
        list(dPhi = pd$dPhi_dlog_thalf, dphi = pd$dphi_dlog_thalf,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dlog_thalf),
        list(dPhi = pd$dPhi_dnu,        dphi = pd$dphi_dnu,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dnu),
        list(dPhi = pd$dPhi_dm,         dphi = pd$dphi_dm,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dm)
      )

      for (s in seq_along(shape_derivs)) {
        dPhi_ds <- shape_derivs[[s]]$dPhi
        dphi_ds <- shape_derivs[[s]]$dphi

        dlogl_ds <- sum(w_H * mu_j * dPhi_ds) +
                    sum(inv_h[idx_event] * mu_j[idx_event] * dphi_ds[idx_event])

        # Left-censoring correction: swap dPhi at time for dPhi at upper
        if (has_left && !is.null(shape_derivs[[s]]$dPhi_upper)) {
          dlogl_ds <- dlogl_ds +
            sum(w_H[idx_left] * mu_j[idx_left] * shape_derivs[[s]]$dPhi_upper) -
            sum(w_H[idx_left] * mu_j[idx_left] * dPhi_ds[idx_left])
        }

        grad[pos] <- dlogl_ds
        pos <- pos + 1L
      }
    } else if (phases[[nm]]$type == "g3") {
      pd <- phase_deriv[[j]]

      # Compute shape derivatives at upper[idx_left] for left-censoring
      pd_upper <- NULL
      if (has_left) {
        tau_j <- exp(pars$log_tau)
        pd_upper <- .hzr_g3_phase_derivatives(upper[idx_left],
          tau = tau_j, gamma = pars$gamma,
          alpha = pars$alpha, eta = pars$eta)
      }

      # For each shape param s in {log_tau, gamma, alpha, eta}:
      shape_derivs <- list(
        list(dPhi = pd$dPhi_dlog_tau, dphi = pd$dphi_dlog_tau,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dlog_tau),
        list(dPhi = pd$dPhi_dgamma,   dphi = pd$dphi_dgamma,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dgamma),
        list(dPhi = pd$dPhi_dalpha,   dphi = pd$dphi_dalpha,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_dalpha),
        list(dPhi = pd$dPhi_deta,     dphi = pd$dphi_deta,
             dPhi_upper = if (!is.null(pd_upper)) pd_upper$dPhi_deta)
      )

      for (s in seq_along(shape_derivs)) {
        dPhi_ds <- shape_derivs[[s]]$dPhi
        dphi_ds <- shape_derivs[[s]]$dphi

        dlogl_ds <- sum(w_H * mu_j * dPhi_ds) +
                    sum(inv_h[idx_event] * mu_j[idx_event] * dphi_ds[idx_event])

        # Left-censoring correction: swap dPhi at time for dPhi at upper
        if (has_left && !is.null(shape_derivs[[s]]$dPhi_upper)) {
          dlogl_ds <- dlogl_ds +
            sum(w_H[idx_left] * mu_j[idx_left] * shape_derivs[[s]]$dPhi_upper) -
            sum(w_H[idx_left] * mu_j[idx_left] * dPhi_ds[idx_left])
        }

        grad[pos] <- dlogl_ds
        pos <- pos + 1L
      }
    }

    # ── Covariate coefficients beta_jk ──────────────────────────────────
    # dH/d(beta_jk) = x_ik * mu_j_i * Phi_j_i
    # dh/d(beta_jk) = x_ik * mu_j_i * phi_j_i
    if (length(pars$beta) > 0 && !is.null(x_list[[nm]])) {
      x_j <- x_list[[nm]]
      for (k in seq_along(pars$beta)) {
        x_k <- x_j[, k]
        dlogl_dbeta <- sum(w_H * x_k * mu_j * Phi_j) +
                       sum(inv_h[idx_event] * x_k[idx_event] *
                           mu_j[idx_event] * phi_j[idx_event])

        if (has_left) {
          dlogl_dbeta <- dlogl_dbeta +
            sum(w_H[idx_left] * x_k[idx_left] *
                mu_j[idx_left] * Phi_j_upper) -
            sum(w_H[idx_left] * x_k[idx_left] *
                mu_j[idx_left] * Phi_j[idx_left])
        }

        grad[pos] <- dlogl_dbeta
        pos <- pos + 1L
      }
    } else {
      pos <- pos + covariate_counts[[nm]]
    }
  }

  # ── Interval-censored fallback: add correction via finite difference ───
  # Interval-censored observations are uncommon; their gradient contribution
  # is added via one-sided numerical difference on the log-likelihood of
  # just those observations.
  if (has_interval) {
    logl_iv <- function(th) {
      cumhaz_l <- .hzr_multiphase_cumhaz(lower, th, phases, covariate_counts,
                                          x_list)
      cumhaz_u <- .hzr_multiphase_cumhaz(upper, th, phases, covariate_counts,
                                          x_list)
      delta <- cumhaz_u[idx_interval] - cumhaz_l[idx_interval]
      sum(-cumhaz_l[idx_interval] + hzr_log1mexp(delta))
    }
    eps_rel <- sqrt(.Machine$double.eps)
    ll0_iv <- logl_iv(theta)
    for (i in seq_len(p)) {
      h_i <- eps_rel * max(abs(theta[i]), 1)
      theta_p <- theta
      theta_p[i] <- theta_p[i] + h_i
      grad[i] <- grad[i] + (logl_iv(theta_p) - ll0_iv) / h_i
    }
  }

  # Safety: zero out non-finite entries
  grad[!is.finite(grad)] <- 0

  grad
}


# ============================================================================
# Optimizer
# ============================================================================

#' Fit a multiphase additive hazard model via maximum likelihood
#'
#' Assembles starting values from phase specifications, resolves per-phase
#' design matrices, and delegates to `.hzr_optim_generic()`.
#'
#' @param time Numeric vector of follow-up times.
#' @param status Numeric event indicator vector.
#' @param time_lower Optional lower bounds for interval censoring.
#' @param time_upper Optional upper bounds for left/interval censoring.
#' @param x Global design matrix (n x p) or NULL.
#' @param theta_start Starting parameter vector (full internal scale).
#'   If NULL, assembled automatically from phase specs.
#' @param control Named list of control options.
#' @param phases Named list of validated `hzr_phase` objects.
#' @param formula_global The global formula (used when phases have no
#'   phase-specific formula).
#' @param data Data frame containing covariates (needed for phase-specific
#'   formula evaluation).
#' @return List with par (internal scale), value, convergence, vcov, etc.
#' @keywords internal
.hzr_optim_multiphase <- function(time, status,
                                   time_lower = NULL, time_upper = NULL,
                                   x = NULL,
                                   theta_start = NULL,
                                   control = list(),
                                   phases,
                                   formula_global = NULL,
                                   data = NULL) {

  phases <- .hzr_validate_phases(phases)

  # --- Resolve per-phase design matrices and covariate counts ----------------
  x_list <- vector("list", length(phases))
  names(x_list) <- names(phases)
  covariate_counts <- setNames(integer(length(phases)), names(phases))

  for (nm in names(phases)) {
    ph <- phases[[nm]]
    if (!is.null(ph$formula) && !is.null(data)) {
      # Phase-specific formula: build design matrix from data
      x_j <- stats::model.matrix(ph$formula, data = data)[, -1L, drop = FALSE]
      x_list[[nm]] <- x_j
      covariate_counts[[nm]] <- ncol(x_j)
    } else if (!is.null(x)) {
      # Inherit global design matrix
      x_list[[nm]] <- x
      covariate_counts[[nm]] <- ncol(x)
    } else {
      x_list[[nm]] <- NULL
      covariate_counts[[nm]] <- 0L
    }
  }

  # --- Assemble starting values if not provided ------------------------------
  if (is.null(theta_start)) {
    theta_start <- unlist(lapply(names(phases), function(nm) {
      .hzr_phase_start(phases[[nm]], n_covariates = covariate_counts[[nm]])
    }))
  }

  # Build parameter names for the full theta vector
  theta_names <- unlist(lapply(names(phases), function(nm) {
    cov_names <- if (covariate_counts[[nm]] > 0 && !is.null(x_list[[nm]])) {
      colnames(x_list[[nm]])
    } else {
      character(0)
    }
    # Ensure we have names even if colnames are NULL
    if (covariate_counts[[nm]] > 0 && length(cov_names) == 0) {
      cov_names <- paste0("x", seq_len(covariate_counts[[nm]]))
    }
    .hzr_phase_theta_names(phases[[nm]], nm, cov_names)
  }))

  names(theta_start) <- theta_names

  # --- Likelihood wrapper matching .hzr_optim_generic() signature -----------
  logl_fn <- function(theta, time, status, time_lower, time_upper, x, ...) {
    .hzr_logl_multiphase(
      theta = theta, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x,
      phases = phases, covariate_counts = covariate_counts,
      x_list = x_list, ...
    )
  }

  # --- Analytic gradient (semi-analytic: chain-rule for mu/beta, central ----
  # differences for shape parameters via .hzr_phase_derivatives()).
  # Falls back to numerical gradient if analytic returns all zeros.
  #
  # IMPORTANT: Save a reference to the ORIGINAL (unwrapped) logl_fn for use

  # inside the gradient function's numerical fallback.  After this point,
  # logl_fn will be re-bound by CoE and fixed-mask wrappers, but the gradient
  # function always receives a FULL theta vector from its own wrapper chain,
  # so it must call the unwrapped version that also expects full theta.
  logl_fn_unwrapped <- logl_fn

  gradient_fn <- function(theta, time, status, time_lower, time_upper, x, ...) {
    grad <- .hzr_gradient_multiphase(
      theta = theta, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper, x = x,
      phases = phases, covariate_counts = covariate_counts, x_list = x_list
    )

    # Fallback: if gradient is all zero (e.g. at infeasible point), try
    # numerical gradient to keep optimizer moving
    if (all(grad == 0)) {
      eps_rel <- sqrt(.Machine$double.eps)
      p <- length(theta)
      ll0 <- logl_fn_unwrapped(theta, time, status, time_lower,
                                time_upper, x, ...)
      for (i in seq_len(p)) {
        h_i <- eps_rel * max(abs(theta[i]), 1)
        theta_plus <- theta
        theta_plus[i] <- theta_plus[i] + h_i
        ll_plus <- logl_fn_unwrapped(theta_plus, time, status,
                                      time_lower, time_upper, x, ...)
        grad[i] <- (ll_plus - ll0) / h_i
      }
    }

    grad
  }

  # --- Conservation of Events (CoE) setup -------------------------------------
  use_conserve <- if (!is.null(control$conserve)) control$conserve else TRUE
  control$conserve <- NULL  # remove before passing to optim

  fixmu_phase <- NULL
  fixmu_pos <- NULL
  total_events <- NULL
  log_mu_positions <- .hzr_log_mu_positions(phases, covariate_counts)

  if (use_conserve && length(phases) >= 2L) {
    total_events <- sum(status == 1)
    if (total_events > 0) {
      # Initial CoE scaling: adjust all log_mu proportionally
      decomp_init <- .hzr_multiphase_cumhaz(
        time, theta_start, phases, covariate_counts, x_list,
        per_phase = TRUE
      )
      sumcz_init <- sum(decomp_init$total)
      if (sumcz_init > 0 && is.finite(sumcz_init)) {
        init_factor <- log(total_events / sumcz_init)
        if (is.finite(init_factor)) {
          for (nm in names(phases)) {
            theta_start[log_mu_positions[[nm]]] <-
              theta_start[log_mu_positions[[nm]]] + init_factor
          }
        }
      }

      # Select fixmu phase (largest cumhaz contributor after scaling)
      fixmu_phase <- .hzr_select_fixmu_phase(
        theta_start, time, status, phases, covariate_counts, x_list
      )
      fixmu_pos <- log_mu_positions[[fixmu_phase]]

      # Apply precise CoE adjustment to the fixmu phase
      theta_start <- .hzr_conserve_events(
        theta_start, fixmu_phase, fixmu_pos,
        time, status, phases, covariate_counts, x_list, total_events
      )

      # Wrap logl and gradient: apply CoE before each evaluation
      logl_fn_pre_coe <- logl_fn
      logl_fn <- function(theta, time, status, time_lower,
                          time_upper, x, ...) {
        theta <- .hzr_conserve_events(
          theta, fixmu_phase, fixmu_pos,
          time, status, phases, covariate_counts, x_list, total_events
        )
        logl_fn_pre_coe(theta, time, status, time_lower,
                        time_upper, x, ...)
      }

      gradient_fn_pre_coe <- gradient_fn
      gradient_fn <- function(theta, time, status, time_lower,
                              time_upper, x, ...) {
        theta <- .hzr_conserve_events(
          theta, fixmu_phase, fixmu_pos,
          time, status, phases, covariate_counts, x_list, total_events
        )
        gradient_fn_pre_coe(theta, time, status, time_lower,
                            time_upper, x, ...)
      }
    } else {
      use_conserve <- FALSE
    }
  } else {
    use_conserve <- FALSE
  }

  # --- Fixed-parameter masking ------------------------------------------------
  # Build a logical mask: TRUE = free (optimized), FALSE = fixed (held constant)
  free_mask <- .hzr_phase_free_mask(phases, covariate_counts)

  # If CoE is active, also fix the fixmu phase's log_mu
  if (use_conserve && !is.null(fixmu_pos)) {
    free_mask[fixmu_pos] <- FALSE
  }

  any_fixed <- !all(free_mask)
  theta_fixed <- theta_start  # full vector; fixed entries stay at these values

  # When some parameters are fixed, wrap logl/gradient to operate on the
  # reduced (free-only) parameter vector.
  if (any_fixed) {
    free_idx <- which(free_mask)

    # Expand reduced theta to full theta
    expand_theta <- function(theta_free) {
      theta_full <- theta_fixed
      theta_full[free_idx] <- theta_free
      theta_full
    }

    logl_fn_full <- logl_fn
    logl_fn <- function(theta, time, status, time_lower, time_upper, x, ...) {
      logl_fn_full(expand_theta(theta), time, status, time_lower,
                   time_upper, x, ...)
    }

    gradient_fn_full <- gradient_fn
    gradient_fn <- function(theta, time, status, time_lower, time_upper, x, ...) {
      grad_full <- gradient_fn_full(expand_theta(theta), time, status,
                                     time_lower, time_upper, x, ...)
      grad_full[free_idx]
    }

    theta_start_optim <- theta_start[free_idx]
  } else {
    theta_start_optim <- theta_start
  }

  # --- Multi-start optimization -----------------------------------------------
  n_starts <- if (!is.null(control$n_starts)) control$n_starts else 5L
  control$n_starts <- NULL  # remove before passing to optim

  best_result <- NULL
  best_value <- -Inf

  # When shapes are fixed the free-parameter count is small (typically just the

  # mu intercepts).  BFGS can struggle on this reduced problem because the
  # likelihood surface is poorly scaled across phases (e.g. the late-phase mu
  # may be ~1e-8).  Strategy: run Nelder-Mead first (derivative-free, robust to
  # scaling) to find a good basin, then polish with BFGS for accurate SE.
  # Nelder-Mead requires >= 2 parameters; use Brent for 1D
  use_nelder_mead_warmup <- any_fixed && length(theta_start_optim) >= 2L &&
    length(theta_start_optim) <= 10L

  for (start_i in seq_len(n_starts)) {
    if (start_i == 1L) {
      theta_try <- theta_start_optim
    } else {
      # Deterministic perturbation seeded from start index for reproducibility
      theta_try <- theta_start_optim +
        stats::rnorm(length(theta_start_optim), sd = 0.5)
    }

    # Nelder-Mead warm-up when shapes are fixed
    if (use_nelder_mead_warmup) {
      # Negative log-likelihood for minimization (matches .hzr_optim_generic)
      nm_objective <- function(theta) {
        if (!all(is.finite(theta))) return(1e10)
        ll <- -logl_fn(
          theta = theta, time = time, status = status,
          time_lower = time_lower, time_upper = time_upper,
          x = x, return_gradient = FALSE
        )
        if (!is.finite(ll)) return(1e10)
        ll
      }
      nm_result <- tryCatch(
        stats::optim(
          par = theta_try, fn = nm_objective,
          method = "Nelder-Mead",
          control = list(maxit = 5000L, reltol = 1e-10)
        ),
        error = function(e) NULL
      )
      if (!is.null(nm_result) && is.finite(nm_result$value)) {
        theta_try <- nm_result$par
      }
    }

    result <- tryCatch(
      .hzr_optim_generic(
        logl_fn     = logl_fn,
        gradient_fn = gradient_fn,
        time        = time, status = status,
        time_lower  = time_lower, time_upper = time_upper,
        x           = x,
        theta_start = theta_try,
        control     = control,
        use_bounds  = FALSE
      ),
      error = function(e) NULL
    )

    if (!is.null(result) && is.finite(result$value) && result$value > best_value) {
      best_value <- result$value
      best_result <- result
    }
  }

  if (is.null(best_result)) {
    stop("Multiphase optimization failed to converge on any start.", call. = FALSE)
  }

  # Expand optimized free params back to full theta vector
  if (any_fixed) {
    theta_full <- theta_fixed
    theta_full[free_idx] <- best_result$par
    best_result$par <- theta_full

    # Final CoE adjustment: solve fixmu log_mu at the optimum
    if (use_conserve && !is.null(fixmu_pos)) {
      best_result$par <- .hzr_conserve_events(
        best_result$par, fixmu_phase, fixmu_pos,
        time, status, phases, covariate_counts, x_list, total_events
      )
    }

    # Expand vcov to full dimension (NA for fixed params — not estimated)
    if (!is.null(best_result$vcov) && is.matrix(best_result$vcov)) {
      p_full <- length(theta_fixed)
      vcov_full <- matrix(NA_real_, p_full, p_full)
      vcov_full[free_idx, free_idx] <- best_result$vcov
      best_result$vcov <- vcov_full
    }

    best_result$fixed_mask <- !free_mask
  }

  # Restore parameter names
  names(best_result$par) <- theta_names

  # Store phase metadata for downstream use (predict, summary)
  best_result$phases <- phases
  best_result$covariate_counts <- covariate_counts
  best_result$x_list <- x_list

  best_result
}


# ============================================================================
# G3 phase derivatives (central finite differences)
# ============================================================================

#' Finite-difference derivatives of G3 phase Phi and phi w.r.t. shape params
#'
#' Computes the G3 cumulative intensity, its time derivative, and their
#' partial derivatives with respect to log_tau, gamma, alpha, and eta using
#' central finite differences.
#'
#' @param time Numeric vector of positive times.
#' @param tau Positive scalar scale parameter.
#' @param gamma Positive scalar time exponent.
#' @param alpha Non-negative scalar shape parameter.
#' @param eta Positive scalar outer exponent.
#' @param h Relative step size for finite differences (default 1e-5).
#' @return Named list with Phi, phi, and 8 derivative vectors.
#' @keywords internal
.hzr_g3_phase_derivatives <- function(time, tau, gamma, alpha, eta,
                                       h = 1e-5) {

  # Base evaluation
  d0 <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                          alpha = alpha, eta = eta)
  Phi0 <- d0$G3
  phi0 <- d0$g3

  # Central differences for log_tau: tau * d/d(tau) = d/d(log(tau))
  eps_tau <- max(abs(tau) * h, 1e-10)
  if (tau - eps_tau > 0) {
    d_plus  <- hzr_decompos_g3(time, tau = tau + eps_tau, gamma = gamma,
                                 alpha = alpha, eta = eta)
    d_minus <- hzr_decompos_g3(time, tau = tau - eps_tau, gamma = gamma,
                                 alpha = alpha, eta = eta)
    dPhi_dtau <- (d_plus$G3 - d_minus$G3) / (2 * eps_tau)
    dphi_dtau <- (d_plus$g3 - d_minus$g3) / (2 * eps_tau)
  } else {
    d_plus  <- hzr_decompos_g3(time, tau = tau + eps_tau, gamma = gamma,
                                 alpha = alpha, eta = eta)
    dPhi_dtau <- (d_plus$G3 - Phi0) / eps_tau
    dphi_dtau <- (d_plus$g3 - phi0) / eps_tau
  }
  # Chain rule: d/d(log_tau) = tau * d/d(tau)
  dPhi_dlog_tau <- tau * dPhi_dtau
  dphi_dlog_tau <- tau * dphi_dtau

  # Central differences for gamma (must stay positive)
  eps_g <- max(abs(gamma) * h, 1e-10)
  if (gamma - eps_g > 0) {
    d_plus  <- hzr_decompos_g3(time, tau = tau, gamma = gamma + eps_g,
                                 alpha = alpha, eta = eta)
    d_minus <- hzr_decompos_g3(time, tau = tau, gamma = gamma - eps_g,
                                 alpha = alpha, eta = eta)
    dPhi_dgamma <- (d_plus$G3 - d_minus$G3) / (2 * eps_g)
    dphi_dgamma <- (d_plus$g3 - d_minus$g3) / (2 * eps_g)
  } else {
    d_plus  <- hzr_decompos_g3(time, tau = tau, gamma = gamma + eps_g,
                                 alpha = alpha, eta = eta)
    dPhi_dgamma <- (d_plus$G3 - Phi0) / eps_g
    dphi_dgamma <- (d_plus$g3 - phi0) / eps_g
  }

  # Central differences for alpha
  if (alpha > h) {
    eps_a <- max(abs(alpha) * h, 1e-10)
    d_plus  <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                 alpha = alpha + eps_a, eta = eta)
    d_minus <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                 alpha = alpha - eps_a, eta = eta)
    dPhi_dalpha <- (d_plus$G3 - d_minus$G3) / (2 * eps_a)
    dphi_dalpha <- (d_plus$g3 - d_minus$g3) / (2 * eps_a)
  } else {
    # alpha near 0: use forward difference
    eps_a <- h
    d_plus <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                alpha = alpha + eps_a, eta = eta)
    dPhi_dalpha <- (d_plus$G3 - Phi0) / eps_a
    dphi_dalpha <- (d_plus$g3 - phi0) / eps_a
  }

  # Central differences for eta (must stay positive)
  eps_e <- max(abs(eta) * h, 1e-10)
  if (eta - eps_e > 0) {
    d_plus  <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                 alpha = alpha, eta = eta + eps_e)
    d_minus <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                 alpha = alpha, eta = eta - eps_e)
    dPhi_deta <- (d_plus$G3 - d_minus$G3) / (2 * eps_e)
    dphi_deta <- (d_plus$g3 - d_minus$g3) / (2 * eps_e)
  } else {
    d_plus  <- hzr_decompos_g3(time, tau = tau, gamma = gamma,
                                 alpha = alpha, eta = eta + eps_e)
    dPhi_deta <- (d_plus$G3 - Phi0) / eps_e
    dphi_deta <- (d_plus$g3 - phi0) / eps_e
  }

  list(
    Phi = Phi0,
    phi = phi0,
    dPhi_dlog_tau = dPhi_dlog_tau,
    dPhi_dgamma   = dPhi_dgamma,
    dPhi_dalpha   = dPhi_dalpha,
    dPhi_deta     = dPhi_deta,
    dphi_dlog_tau = dphi_dlog_tau,
    dphi_dgamma   = dphi_dgamma,
    dphi_dalpha   = dphi_dalpha,
    dphi_deta     = dphi_deta
  )
}


# ============================================================================
# Back-transformation helpers
# ============================================================================

#' Transform multiphase theta from internal to user-facing scale
#'
#' Converts log_mu -> mu, log_t_half -> t_half; nu, m, beta pass through.
#'
#' @param theta Full parameter vector (internal scale).
#' @param phases Named list of validated `hzr_phase` objects.
#' @param covariate_counts Named integer vector of per-phase covariate counts.
#' @return Named numeric vector on the reporting scale.
#' @keywords internal
.hzr_multiphase_theta_user <- function(theta, phases, covariate_counts) {
  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  values <- numeric(0)
  nms    <- character(0)

  for (nm in names(phases)) {
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])

    if (phases[[nm]]$type == "constant") {
      values <- c(values, exp(pars$log_mu))
      nms    <- c(nms, paste0(nm, ".mu"))
    } else if (phases[[nm]]$type == "g3") {
      values <- c(values, exp(pars$log_mu), exp(pars$log_tau),
                  pars$gamma, pars$alpha, pars$eta)
      nms    <- c(nms, paste0(nm, ".mu"), paste0(nm, ".tau"),
                  paste0(nm, ".gamma"), paste0(nm, ".alpha"), paste0(nm, ".eta"))
    } else {
      values <- c(values, exp(pars$log_mu), exp(pars$log_t_half),
                  pars$nu, pars$m)
      nms    <- c(nms, paste0(nm, ".mu"), paste0(nm, ".t_half"),
                  paste0(nm, ".nu"), paste0(nm, ".m"))
    }

    if (length(pars$beta) > 0) {
      values <- c(values, pars$beta)
      # Try to derive covariate names from the theta_split names
      n_skip <- switch(phases[[nm]]$type, constant = 1L, g3 = 5L, 4L)
      raw_names <- names(theta_split[[nm]])
      if (!is.null(raw_names) && length(raw_names) > n_skip) {
        beta_labels <- raw_names[(n_skip + 1L):length(raw_names)]
        # Strip phase prefix if present
        beta_labels <- sub(paste0("^", nm, "\\."), "", beta_labels)
      } else {
        beta_labels <- paste0("x", seq_along(pars$beta))
      }
      nms <- c(nms, paste0(nm, ".", beta_labels))
    }
  }

  names(values) <- nms
  values
}
