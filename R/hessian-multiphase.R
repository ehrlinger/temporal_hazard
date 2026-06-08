#' @keywords internal
NULL

# hessian-multiphase.R -- Analytic NLL Hessian for multiphase hazard models.
#
# Three helpers then the main function:
#   .hzr_phase_second_derivatives()     -- d2Phi/dss' and d2phi/dss' for CDF/hazard
#   .hzr_g3_phase_second_derivatives()  -- same for G3 (late) phases
#   .hzr_hessian_multiphase()           -- full p x p NLL Hessian

# Step size for 4-point second-difference formulas: h ~ eps^(1/4).
# Optimal for central second differences.
.hzr_h2 <- .Machine$double.eps^(1 / 4)

# ---------------------------------------------------------------------------
# CDF / hazard phase second derivatives
# ---------------------------------------------------------------------------

#' Second derivatives of CDF/hazard phase shape functions
#'
#' Computes the six unique mixed second partial derivatives of
#' \eqn{\Phi_j(t)} (cumulative shape) and \eqn{\phi_j(t)} (instantaneous
#' shape) with respect to the internal shape parameters
#' \eqn{\{log\_t\_half, \nu, m\}} via second-order central differences.
#'
#' This is the second-derivative analogue of \code{.hzr_phase_derivatives()}.
#'
#' @param time  Numeric vector of observation times.
#' @param t_half,nu,m  Phase shape parameters (natural scale).
#' @param type  Phase type: \code{"cdf"} or \code{"hazard"}.
#' @return A list with elements \code{d2Phi_dlog_thalf2},
#'   \code{d2Phi_dnu2}, \code{d2Phi_dm2},
#'   \code{d2Phi_dlog_thalf_dnu}, \code{d2Phi_dlog_thalf_dm},
#'   \code{d2Phi_dnu_dm}, and the matching \code{d2phi_*} entries.
#' @noRd
.hzr_phase_second_derivatives <- function(time, t_half, nu, m, type) {
  h <- .hzr_h2

  # Evaluate Phi and phi at given (t_half_p, nu_p, m_p)
  eval_at <- function(th, nup, mp) {
    d <- .hzr_phase_derivatives(time, t_half = th, nu = nup, m = mp, type = type)
    list(Phi = d$Phi, phi = d$phi)
  }

  d00 <- eval_at(t_half, nu, m)

  # Perturbed values on each internal scale:
  #   log_t_half: steps are multiplicative -> t_half * exp(±h)
  #   nu, m:      additive steps
  t_p <- t_half * exp( h);  t_m <- t_half * exp(-h)
  nu_p <- nu + h;            nu_m <- nu - h
  m_p  <- m  + h;            m_m  <- m  - h

  # Six single-parameter perturbations for diagonal second differences
  d_tp  <- eval_at(t_p,    nu,   m)
  d_tm  <- eval_at(t_m,    nu,   m)
  d_np  <- eval_at(t_half, nu_p, m)
  d_nm  <- eval_at(t_half, nu_m, m)
  d_mp  <- eval_at(t_half, nu,   m_p)
  d_mm  <- eval_at(t_half, nu,   m_m)

  h2 <- h * h

  # Diagonal second differences
  d2Phi_tt  <- (d_tp$Phi - 2 * d00$Phi + d_tm$Phi) / h2
  d2Phi_nn  <- (d_np$Phi - 2 * d00$Phi + d_nm$Phi) / h2
  d2Phi_mm_ <- (d_mp$Phi - 2 * d00$Phi + d_mm$Phi) / h2
  d2phi_tt  <- (d_tp$phi - 2 * d00$phi + d_tm$phi) / h2
  d2phi_nn  <- (d_np$phi - 2 * d00$phi + d_nm$phi) / h2
  d2phi_mm_ <- (d_mp$phi - 2 * d00$phi + d_mm$phi) / h2

  # Four-point mixed second differences
  # d2f/dada db = [f(a+h,b+h) - f(a+h,b-h) - f(a-h,b+h) + f(a-h,b-h)] / (4h^2)
  four_pt_Phi <- function(p1p, p1m, p2p, p2m) {
    dpp <- eval_at(p1p[[1]], p1p[[2]], p1p[[3]])
    dpm <- eval_at(p1p[[1]], p1p[[2]], p1p[[3]])  # placeholder, overridden
    dmp <- eval_at(p1m[[1]], p1m[[2]], p1m[[3]])
    dmm <- eval_at(p1m[[1]], p1m[[2]], p1m[[3]])  # placeholder, overridden
    # Need all four corners: (p1+, p2+), (p1+, p2-), (p1-, p2+), (p1-, p2-)
    e_pp <- eval_at(p1p[[1]], p2p[[2]], p2p[[3]])
    e_pm <- eval_at(p1p[[1]], p2m[[2]], p2m[[3]])
    e_mp <- eval_at(p1m[[1]], p2p[[2]], p2p[[3]])
    e_mm <- eval_at(p1m[[1]], p2m[[2]], p2m[[3]])
    list(
      Phi = (e_pp$Phi - e_pm$Phi - e_mp$Phi + e_mm$Phi) / (4 * h2),
      phi = (e_pp$phi - e_pm$phi - e_mp$phi + e_mm$phi) / (4 * h2)
    )
  }

  # t_half / nu cross: corners (t_p, nu_p), (t_p, nu_m), (t_m, nu_p), (t_m, nu_m)
  e_tn_pp <- eval_at(t_p, nu_p, m);  e_tn_pm <- eval_at(t_p, nu_m, m)
  e_tn_mp <- eval_at(t_m, nu_p, m);  e_tn_mm <- eval_at(t_m, nu_m, m)
  d2_tn_Phi <- (e_tn_pp$Phi - e_tn_pm$Phi - e_tn_mp$Phi + e_tn_mm$Phi) / (4 * h2)
  d2_tn_phi <- (e_tn_pp$phi - e_tn_pm$phi - e_tn_mp$phi + e_tn_mm$phi) / (4 * h2)

  # t_half / m cross
  e_tm_pp <- eval_at(t_p, nu, m_p);  e_tm_pm <- eval_at(t_p, nu, m_m)
  e_tm_mp <- eval_at(t_m, nu, m_p);  e_tm_mm <- eval_at(t_m, nu, m_m)
  d2_tm_Phi <- (e_tm_pp$Phi - e_tm_pm$Phi - e_tm_mp$Phi + e_tm_mm$Phi) / (4 * h2)
  d2_tm_phi <- (e_tm_pp$phi - e_tm_pm$phi - e_tm_mp$phi + e_tm_mm$phi) / (4 * h2)

  # nu / m cross
  e_nm_pp <- eval_at(t_half, nu_p, m_p);  e_nm_pm <- eval_at(t_half, nu_p, m_m)
  e_nm_mp <- eval_at(t_half, nu_m, m_p);  e_nm_mm <- eval_at(t_half, nu_m, m_m)
  d2_nm_Phi <- (e_nm_pp$Phi - e_nm_pm$Phi - e_nm_mp$Phi + e_nm_mm$Phi) / (4 * h2)
  d2_nm_phi <- (e_nm_pp$phi - e_nm_pm$phi - e_nm_mp$phi + e_nm_mm$phi) / (4 * h2)

  list(
    d2Phi_dlog_thalf2    = d2Phi_tt,
    d2Phi_dnu2           = d2Phi_nn,
    d2Phi_dm2            = d2Phi_mm_,
    d2Phi_dlog_thalf_dnu = d2_tn_Phi,
    d2Phi_dlog_thalf_dm  = d2_tm_Phi,
    d2Phi_dnu_dm         = d2_nm_Phi,

    d2phi_dlog_thalf2    = d2phi_tt,
    d2phi_dnu2           = d2phi_nn,
    d2phi_dm2            = d2phi_mm_,
    d2phi_dlog_thalf_dnu = d2_tn_phi,
    d2phi_dlog_thalf_dm  = d2_tm_phi,
    d2phi_dnu_dm         = d2_nm_phi
  )
}

# ---------------------------------------------------------------------------
# G3 (late) phase second derivatives
# ---------------------------------------------------------------------------

#' Second derivatives of G3 phase shape functions
#'
#' Analogue of \code{.hzr_g3_phase_derivatives()} for second-order
#' shape derivatives over \code{\{log\_tau, gamma, alpha, eta\}}.
#'
#' @param time,tau,gamma,alpha,eta  G3 phase parameters.
#' @return List with \code{d2Phi_*} and \code{d2phi_*} entries for the ten
#'   unique pairs from \code{\{log\_tau, gamma, alpha, eta\}}.
#' @noRd
.hzr_g3_phase_second_derivatives <- function(time, tau, gamma, alpha, eta) {
  h <- .hzr_h2

  eval_at <- function(tau_p, gam_p, alp_p, eta_p) {
    # Guard against non-positive G3 parameters
    if (tau_p <= 0 || gam_p <= 0 || eta_p <= 0) return(NULL)
    tryCatch(
      {
        d <- hzr_decompos_g3(time, tau = tau_p, gamma = gam_p,
                              alpha = alp_p, eta = eta_p)
        list(Phi = d$G3, phi = d$g3)
      },
      error = function(e) NULL
    )
  }

  d00 <- eval_at(tau, gamma, alpha, eta)
  if (is.null(d00)) {
    n <- length(time)
    zeros <- rep(0, n)
    nms <- c("d2Phi_dlog_tau2", "d2Phi_dgamma2", "d2Phi_dalpha2",
             "d2Phi_deta2",
             "d2Phi_dlog_tau_dgamma", "d2Phi_dlog_tau_dalpha",
             "d2Phi_dlog_tau_deta",   "d2Phi_dgamma_dalpha",
             "d2Phi_dgamma_deta",     "d2Phi_dalpha_deta",
             "d2phi_dlog_tau2", "d2phi_dgamma2", "d2phi_dalpha2",
             "d2phi_deta2",
             "d2phi_dlog_tau_dgamma", "d2phi_dlog_tau_dalpha",
             "d2phi_dlog_tau_deta",   "d2phi_dgamma_dalpha",
             "d2phi_dgamma_deta",     "d2phi_dalpha_deta")
    return(setNames(replicate(length(nms), zeros, simplify = FALSE), nms))
  }

  # Perturbed values: log_tau, gamma, alpha, eta steps
  tau_p <- tau * exp( h);  tau_m <- tau * exp(-h)
  gam_p <- gamma + h;      gam_m <- gamma - h
  alp_p <- alpha + h;      alp_m <- alpha - h
  eta_p <- eta   + h;      eta_m <- eta   - h

  # Fallback to d00 if perturbed eval returns NULL
  safe_eval <- function(tau_v, gam_v, alp_v, eta_v) {
    r <- eval_at(tau_v, gam_v, alp_v, eta_v)
    if (is.null(r)) d00 else r
  }

  d_Tp <- safe_eval(tau_p, gamma, alpha, eta)
  d_Tm <- safe_eval(tau_m, gamma, alpha, eta)
  d_Gp <- safe_eval(tau, gam_p, alpha, eta)
  d_Gm <- safe_eval(tau, gam_m, alpha, eta)
  d_Ap <- safe_eval(tau, gamma, alp_p, eta)
  d_Am <- safe_eval(tau, gamma, alp_m, eta)
  d_Ep <- safe_eval(tau, gamma, alpha, eta_p)
  d_Em <- safe_eval(tau, gamma, alpha, eta_m)

  h2 <- h * h

  # Diagonal second differences
  diag2 <- function(dp, dm, d0, which) {
    (dp[[which]] - 2 * d0[[which]] + dm[[which]]) / h2
  }

  # Four-point mixed helper
  mixed2 <- function(tau_a, gam_a, alp_a, eta_a,
                     tau_b, gam_b, alp_b, eta_b, which) {
    e_pp <- safe_eval(tau_a, gam_a, alp_a, eta_a)
    e_pm <- safe_eval(tau_a, gam_a, alp_a, eta_a)  # overridden below
    e_mp <- safe_eval(tau_b, gam_b, alp_b, eta_b)  # overridden below
    e_mm <- safe_eval(tau_b, gam_b, alp_b, eta_b)  # overridden below
    # Correct: we need (par1+, par2+), (par1+, par2-), (par1-, par2+), (par1-, par2-)
    # Since the two parameters are independent, par1+/par2+ is the combination
    # This is handled directly below for each pair
    stop("use explicit four-point calls")
  }

  # tau / gamma
  d_TpGp <- safe_eval(tau_p, gam_p, alpha, eta)
  d_TpGm <- safe_eval(tau_p, gam_m, alpha, eta)
  d_TmGp <- safe_eval(tau_m, gam_p, alpha, eta)
  d_TmGm <- safe_eval(tau_m, gam_m, alpha, eta)
  TG_Phi <- (d_TpGp$Phi - d_TpGm$Phi - d_TmGp$Phi + d_TmGm$Phi) / (4 * h2)
  TG_phi <- (d_TpGp$phi - d_TpGm$phi - d_TmGp$phi + d_TmGm$phi) / (4 * h2)

  # tau / alpha
  d_TpAp <- safe_eval(tau_p, gamma, alp_p, eta)
  d_TpAm <- safe_eval(tau_p, gamma, alp_m, eta)
  d_TmAp <- safe_eval(tau_m, gamma, alp_p, eta)
  d_TmAm <- safe_eval(tau_m, gamma, alp_m, eta)
  TA_Phi <- (d_TpAp$Phi - d_TpAm$Phi - d_TmAp$Phi + d_TmAm$Phi) / (4 * h2)
  TA_phi <- (d_TpAp$phi - d_TpAm$phi - d_TmAp$phi + d_TmAm$phi) / (4 * h2)

  # tau / eta
  d_TpEp <- safe_eval(tau_p, gamma, alpha, eta_p)
  d_TpEm <- safe_eval(tau_p, gamma, alpha, eta_m)
  d_TmEp <- safe_eval(tau_m, gamma, alpha, eta_p)
  d_TmEm <- safe_eval(tau_m, gamma, alpha, eta_m)
  TE_Phi <- (d_TpEp$Phi - d_TpEm$Phi - d_TmEp$Phi + d_TmEm$Phi) / (4 * h2)
  TE_phi <- (d_TpEp$phi - d_TpEm$phi - d_TmEp$phi + d_TmEm$phi) / (4 * h2)

  # gamma / alpha
  d_GpAp <- safe_eval(tau, gam_p, alp_p, eta)
  d_GpAm <- safe_eval(tau, gam_p, alp_m, eta)
  d_GmAp <- safe_eval(tau, gam_m, alp_p, eta)
  d_GmAm <- safe_eval(tau, gam_m, alp_m, eta)
  GA_Phi <- (d_GpAp$Phi - d_GpAm$Phi - d_GmAp$Phi + d_GmAm$Phi) / (4 * h2)
  GA_phi <- (d_GpAp$phi - d_GpAm$phi - d_GmAp$phi + d_GmAm$phi) / (4 * h2)

  # gamma / eta
  d_GpEp <- safe_eval(tau, gam_p, alpha, eta_p)
  d_GpEm <- safe_eval(tau, gam_p, alpha, eta_m)
  d_GmEp <- safe_eval(tau, gam_m, alpha, eta_p)
  d_GmEm <- safe_eval(tau, gam_m, alpha, eta_m)
  GE_Phi <- (d_GpEp$Phi - d_GpEm$Phi - d_GmEp$Phi + d_GmEm$Phi) / (4 * h2)
  GE_phi <- (d_GpEp$phi - d_GpEm$phi - d_GmEp$phi + d_GmEm$phi) / (4 * h2)

  # alpha / eta
  d_ApEp <- safe_eval(tau, gamma, alp_p, eta_p)
  d_ApEm <- safe_eval(tau, gamma, alp_p, eta_m)
  d_AmEp <- safe_eval(tau, gamma, alp_m, eta_p)
  d_AmEm <- safe_eval(tau, gamma, alp_m, eta_m)
  AE_Phi <- (d_ApEp$Phi - d_ApEm$Phi - d_AmEp$Phi + d_AmEm$Phi) / (4 * h2)
  AE_phi <- (d_ApEp$phi - d_ApEm$phi - d_AmEp$phi + d_AmEm$phi) / (4 * h2)

  list(
    d2Phi_dlog_tau2        = diag2(d_Tp, d_Tm, d00, "Phi"),
    d2Phi_dgamma2          = diag2(d_Gp, d_Gm, d00, "Phi"),
    d2Phi_dalpha2          = diag2(d_Ap, d_Am, d00, "Phi"),
    d2Phi_deta2            = diag2(d_Ep, d_Em, d00, "Phi"),
    d2Phi_dlog_tau_dgamma  = TG_Phi,
    d2Phi_dlog_tau_dalpha  = TA_Phi,
    d2Phi_dlog_tau_deta    = TE_Phi,
    d2Phi_dgamma_dalpha    = GA_Phi,
    d2Phi_dgamma_deta      = GE_Phi,
    d2Phi_dalpha_deta      = AE_Phi,

    d2phi_dlog_tau2        = diag2(d_Tp, d_Tm, d00, "phi"),
    d2phi_dgamma2          = diag2(d_Gp, d_Gm, d00, "phi"),
    d2phi_dalpha2          = diag2(d_Ap, d_Am, d00, "phi"),
    d2phi_deta2            = diag2(d_Ep, d_Em, d00, "phi"),
    d2phi_dlog_tau_dgamma  = TG_phi,
    d2phi_dlog_tau_dalpha  = TA_phi,
    d2phi_dlog_tau_deta    = TE_phi,
    d2phi_dgamma_dalpha    = GA_phi,
    d2phi_dgamma_deta      = GE_phi,
    d2phi_dalpha_deta      = AE_phi
  )
}

# ---------------------------------------------------------------------------
# Main analytic Hessian function (.hzr_hessian_multiphase)
# ---------------------------------------------------------------------------
# (Added in Task 2 below)
