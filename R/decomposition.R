# decomposition.R — Generalized temporal decomposition
#
# PURPOSE
# -------
# The unified parametric family decompos(t; t_half, nu, m) generates all
# temporal phase shapes used in multiphase hazard models.  It produces three
# quantities: G(t) (CDF), g(t) (density — "early" pattern), and h(t) =
# g(t)/(1-G(t)) (hazard — "late" pattern).
#
# Originally introduced by Blackstone, Naftel, and Turner (1986, JASA 81:615)
# and extended to longitudinal mixed-effects settings by Rajeswaran et al.
# (2018, Stat Methods Med Res 27:126).  Ported here from the mixhazard
# package with enhanced numerical stability and integrated phase-level
# helpers for the additive cumulative hazard model:
#
#   H(t | x) = sum_j  mu_j(x) * Phi_j(t; t_half_j, nu_j, m_j)
#
# PARAMETER MAPPING FROM SAS/C HAZARD
# ------------------------------------
# SAS Early (G1): DELTA, RHO/THALF, NU, M  -->  t_half, nu, m
# SAS Late  (G3): TAU, GAMMA, ALPHA, ETA   -->  t_half, nu, m
# SAS Const (G2): (none)                   -->  type = "constant"
#
# The 4-parameter C/SAS parameterizations collapse onto this 3-parameter
# family.  The C DELTA parameter controlled a time transformation
# B(t) = (exp(delta*t) - 1)/delta that is absorbed by the shape.

# ============================================================================
# hzr_decompos — core engine
# ============================================================================

#' Generalized temporal decomposition
#'
#' Computes the cumulative distribution \eqn{G(t)}, density \eqn{g(t)}, and
#' hazard \eqn{h(t) = g(t)/(1 - G(t))} for the parametric family defined by
#' half-life, time exponent, and shape.  This single function generates all
#' temporal phase shapes used in multiphase hazard models.
#'
#' @section Parameter mapping from SAS/C HAZARD:
#'
#' The original C code used separate parameterizations for early (DELTA,
#' RHO/THALF, NU, M) and late (TAU, GAMMA, ALPHA, ETA) phases.  Both
#' collapse onto the three parameters here.  See
#' [hzr_argument_mapping()] for the full translation table.
#'
#' @section Valid parameter combinations:
#'
#' Six cases are defined by the signs of `nu` and `m`:
#'
#' \tabular{lll}{
#'   **Case** \tab **Sign** \tab **Behavior** \cr
#'   1     \tab m > 0, nu > 0 \tab Standard sigmoidal \cr
#'   1L    \tab m = 0, nu > 0 \tab Exponential-like (Weibull CDF) \cr
#'   2     \tab m < 0, nu > 0 \tab Heavy-tailed \cr
#'   2L    \tab m < 0, nu = 0 \tab Exponential decay \cr
#'   3     \tab m > 0, nu < 0 \tab Bounded cumulative \cr
#'   3L    \tab m = 0, nu < 0 \tab Bounded exponential \cr
#' }
#'
#' The combination m < 0 **and** nu < 0 is undefined and raises an error.
#'
#' @param time Numeric vector of times (must be > 0).
#' @param t_half Half-life: time at which \eqn{G(t_{1/2}) = 0.5}.
#'   Must be > 0.
#' @param nu Time exponent controlling rate dynamics.
#'   SAS early: `NU`.  SAS late: relates to `GAMMA`/`ETA`.
#' @param m Shape exponent controlling the distributional form.
#'   SAS early: `M`.  SAS late: relates to `GAMMA`/`ALPHA`.
#'
#' @return A named list with three numeric vectors, each the same length
#'   as `time`:
#' \describe{
#'   \item{G}{Cumulative distribution \eqn{G(t) \in [0, 1]}.}
#'   \item{g}{Density \eqn{g(t) = dG/dt \ge 0}.  The "early" phase
#'     temporal pattern.}
#'   \item{h}{Hazard \eqn{h(t) = g(t)/(1 - G(t)) \ge 0}.  The "late"
#'     phase temporal pattern.}
#' }
#'
#' @references
#' Blackstone EH, Naftel DC, Turner ME Jr. The decomposition of time-varying
#' hazard into phases, each incorporating a separate stream of concomitant
#' information. *J Am Stat Assoc.* 1986;81(395):615--624.
#' \doi{10.1080/01621459.1986.10478314}
#'
#' Rajeswaran J, Blackstone EH, Ehrlinger J, Li L, Ishwaran H, Parides MK.
#' Probability of atrial fibrillation after ablation: Using a parametric
#' nonlinear temporal decomposition mixed effects model. *Stat Methods Med Res.*
#' 2018;27(1):126--141. \doi{10.1177/0962280215623583}
#'
#' @examples
#' t_grid <- seq(0.1, 10, by = 0.1)
#'
#' # Case 1: standard sigmoidal (m > 0, nu > 0)
#' d1 <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 1)
#' plot(t_grid, d1$G, type = "l", main = "CDF (m=1, nu=2)")
#'
#' # Case 1L: Weibull-like (m = 0, nu > 0)
#' d1L <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
#'
#' # Case 2: heavy-tailed (m < 0, nu > 0)
#' d2 <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = -1)
#'
#' @seealso [hzr_phase_cumhaz()] for the phase-level cumulative hazard
#'   contribution, [hzr_argument_mapping()] for SAS/C parameter mapping.
#'
#' @export
hzr_decompos <- function(time, t_half, nu, m) {

  # --- Input validation ------------------------------------------------------
  stopifnot(
    "time must be numeric"   = is.numeric(time),
    "t_half must be a positive scalar" =
      is.numeric(t_half) && length(t_half) == 1L && t_half > 0,
    "nu must be a numeric scalar" =
      is.numeric(nu) && length(nu) == 1L && is.finite(nu),
    "m must be a numeric scalar" =
      is.numeric(m) && length(m) == 1L && is.finite(m)
  )

  if (m < 0 && nu < 0) {
    stop("Decomposition undefined when both m < 0 and nu < 0 (m = ",
         m, ", nu = ", nu, ").", call. = FALSE)
  }

  # Clamp time away from zero to avoid division-by-zero / 0^negative

  time <- pmax(time, .Machine$double.xmin)

  n <- length(time)

  # Pre-compute recurring exponent terms
  if (m != 0) mm1 <- -(1 / m) - 1
  if (nu != 0) num1 <- -(1 / nu) - 1

  # --- Case dispatch ---------------------------------------------------------
  if (m > 0 && nu > 0) {
    # Case 1: standard sigmoidal
    rho   <- nu * t_half * (((2^m - 1) / m)^nu)
    bt    <- nu * time / rho
    btnu  <- 1 + m * bt^(-1 / nu)
    G     <- btnu^(-1 / m)
    g     <- (btnu^mm1) * (bt^num1) / rho

  } else if (m == 0 && nu > 0) {
    # Case 1L: Weibull-like (m -> 0 limit)
    rho   <- nu * t_half * (log(2)^nu)
    bt    <- nu * time / rho
    btnu  <- bt^(-1 / nu)
    G     <- exp(-btnu)
    g     <- G * (bt^num1) / rho

  } else if (m < 0 && nu > 0) {
    # Case 2: heavy-tailed
    rho   <- nu * t_half / ((1 - 2^m)^(-nu) - 1)
    bt    <- 1 + nu * time / rho
    btnu  <- 1 - bt^(-1 / nu)
    G     <- btnu^(-1 / m)
    g     <- -(btnu^mm1) * (bt^num1) / (m * rho)

  } else if (m < 0 && nu == 0) {
    # Case 2L: exponential decay (nu -> 0 limit)
    rho   <- -t_half / log(1 - 2^m)
    bt    <- exp(-time / rho)
    btm   <- 1 - bt
    G     <- btm^(-1 / m)
    g     <- -(btm^mm1) * bt / (m * rho)

  } else if (m > 0 && nu < 0) {
    # Case 3: bounded cumulative
    rho   <- -nu * t_half * ((2^m - 1)^nu)
    bt    <- -nu * time / rho
    btnu  <- 1 + m * bt^(-1 / nu)
    G     <- 1 - btnu^(-1 / m)
    g     <- (btnu^mm1) * (bt^num1) / rho

  } else if (m == 0 && nu < 0) {
    # Case 3L: bounded exponential (m -> 0 limit)
    rho   <- -nu * t_half * (log(2)^nu)
    bt    <- -nu * time / rho
    btnu  <- bt^(-1 / nu)
    G     <- 1 - exp(-btnu)
    g     <- exp(-btnu) * (bt^num1) / rho
  }

  # --- Hazard from density and CDF ------------------------------------------
  # h(t) = g(t) / (1 - G(t)), with guard against G(t) = 1
  one_minus_G <- pmax(1 - G, .Machine$double.xmin)
  h <- g / one_minus_G

  list(G = G, g = g, h = h)
}


# ============================================================================
# Phase-level helpers for the additive cumulative hazard model
# ============================================================================

#' Cumulative hazard contribution from a single phase
#'
#' Computes \eqn{\Phi_j(t)} for one phase in the additive model
#' \eqn{H(t|x) = \sum_j \mu_j(x) \Phi_j(t)}.
#'
#' @param time Numeric vector of times (> 0).
#' @param t_half Half-life parameter (> 0).
#' @param nu Time exponent.
#' @param m Shape parameter.
#' @param type Phase type: `"cdf"` (early — uses \eqn{G(t)}),
#'   `"hazard"` (late — uses cumulative hazard from \eqn{h(t)}), or
#'   `"constant"` (flat rate — \eqn{\Phi = t}).
#'
#' @return Numeric vector of cumulative hazard contributions \eqn{\Phi(t)},
#'   same length as `time`.
#'
#' @details
#' - `"cdf"`: \eqn{\Phi(t) = G(t)}.  Bounded \eqn{[0, 1]}.  Models early
#'   risk that resolves over time.
#' - `"hazard"`: \eqn{\Phi(t) = -\log(1 - G(t))}.  Monotone increasing.
#'   Models late or aging risk.  This is the cumulative hazard derived from
#'   the hazard function \eqn{h(t)}, since
#'   \eqn{\int_0^t h(s)\,ds = -\log(1 - G(t))}.
#' - `"constant"`: \eqn{\Phi(t) = t}.  Ignores `t_half`, `nu`, `m`.
#'   Equivalent to exponential (constant hazard rate).
#'
#' @examples
#' t_grid <- seq(0.1, 10, by = 0.1)
#' phi_early <- hzr_phase_cumhaz(t_grid, t_half = 2, nu = 2, m = 0,
#'                                type = "cdf")
#' phi_late  <- hzr_phase_cumhaz(t_grid, t_half = 5, nu = 1, m = 0,
#'                                type = "hazard")
#' phi_const <- hzr_phase_cumhaz(t_grid, type = "constant")
#'
#' @seealso [hzr_decompos()] for the underlying parametric family,
#'   [hzr_phase_hazard()] for the instantaneous hazard contribution.
#'
#' @export
hzr_phase_cumhaz <- function(time, t_half = 1, nu = 1, m = 0,
                              type = c("cdf", "hazard", "constant")) {
  type <- match.arg(type)

  if (type == "constant") {
    return(time)
  }

  d <- hzr_decompos(time, t_half = t_half, nu = nu, m = m)

  switch(type,
    cdf    = d$G,
    hazard = {
      # Cumulative hazard from h(t): integral_0^t h(s)ds = -log(1 - G(t))
      one_minus_G <- pmax(1 - d$G, .Machine$double.xmin)
      -log(one_minus_G)
    }
  )
}


# ============================================================================
# Phase-level derivatives for analytic gradient
# ============================================================================

#' Derivatives of phase cumulative and instantaneous hazard w.r.t. shape params
#'
#' Computes \eqn{\Phi_j(t)}, \eqn{\phi_j(t)}, and their derivatives with
#' respect to `t_half`, `nu`, and `m` using central finite differences on
#' [hzr_decompos()].  The `log_t_half` derivative is obtained via the chain
#' rule: \eqn{d\Phi/d(\log t_{1/2}) = t_{1/2} \cdot d\Phi/dt_{1/2}}.
#'
#' @param time Numeric vector of positive times.
#' @param t_half Positive scalar half-life.
#' @param nu Numeric scalar time exponent.
#' @param m Numeric scalar shape exponent.
#' @param type Phase type: `"cdf"`, `"hazard"`, or `"constant"`.
#'
#' @return Named list:
#' \describe{
#'   \item{Phi}{Cumulative hazard contribution \eqn{\Phi(t)}.}
#'   \item{phi}{Instantaneous hazard contribution \eqn{\phi(t) = d\Phi/dt}.}
#'   \item{dPhi_dlog_thalf}{\eqn{d\Phi / d(\log t_{1/2})}.}
#'   \item{dPhi_dnu}{\eqn{d\Phi / d\nu}.}
#'   \item{dPhi_dm}{\eqn{d\Phi / dm}.}
#'   \item{dphi_dlog_thalf}{\eqn{d\phi / d(\log t_{1/2})}.}
#'   \item{dphi_dnu}{\eqn{d\phi / d\nu}.}
#'   \item{dphi_dm}{\eqn{d\phi / dm}.}
#' }
#' @keywords internal
.hzr_phase_derivatives <- function(time, t_half, nu, m,
                                    type = c("cdf", "hazard", "constant")) {
  type <- match.arg(type)
  n <- length(time)

  # Constant phases: Phi = t, phi = 1, no shape dependence

  if (type == "constant") {
    zeros <- rep(0, n)
    return(list(
      Phi = time,
      phi = rep(1, n),
      dPhi_dlog_thalf = zeros,
      dPhi_dnu        = zeros,
      dPhi_dm         = zeros,
      dphi_dlog_thalf = zeros,
      dphi_dnu        = zeros,
      dphi_dm         = zeros
    ))
  }

  # Base evaluation
  d0 <- hzr_decompos(time, t_half = t_half, nu = nu, m = m)

  # Extract Phi and phi from decomposition for the given phase type
  extract <- function(d, tp) {
    if (tp == "cdf") {
      Phi <- d$G
      phi <- d$g
    } else {
      # "hazard": Phi = -log(1 - G), phi = h
      one_minus_G <- pmax(1 - d$G, .Machine$double.xmin)
      Phi <- -log(one_minus_G)
      phi <- d$h
    }
    list(Phi = Phi, phi = phi)
  }

  base <- extract(d0, type)

  # Central-difference derivatives w.r.t. shape parameters
  eps_rel <- (.Machine$double.eps)^(1/3)  # optimal for central diff

  # Helper: perturbed decomposition (returns NULL if invalid params)
  perturb_decompos <- function(th, n_val, m_val) {
    if (m_val < 0 && n_val < 0) return(NULL)
    if (th <= 0) return(NULL)
    tryCatch(
      hzr_decompos(time, t_half = th, nu = n_val, m = m_val),
      error = function(e) NULL
    )
  }

  # Derivative w.r.t. t_half (then multiply by t_half for log_t_half chain rule)
  h_th <- eps_rel * max(abs(t_half), 1e-4)
  d_plus  <- perturb_decompos(t_half + h_th, nu, m)
  d_minus <- perturb_decompos(t_half - h_th, nu, m)
  if (!is.null(d_plus) && !is.null(d_minus)) {
    e_plus  <- extract(d_plus, type)
    e_minus <- extract(d_minus, type)
    dPhi_dt_half <- (e_plus$Phi - e_minus$Phi) / (2 * h_th)
    dphi_dt_half <- (e_plus$phi - e_minus$phi) / (2 * h_th)
  } else {
    # One-sided fallback
    if (!is.null(d_plus)) {
      e_plus <- extract(d_plus, type)
      dPhi_dt_half <- (e_plus$Phi - base$Phi) / h_th
      dphi_dt_half <- (e_plus$phi - base$phi) / h_th
    } else {
      dPhi_dt_half <- rep(0, n)
      dphi_dt_half <- rep(0, n)
    }
  }
  # Chain rule: d/d(log_t_half) = t_half * d/d(t_half)
  dPhi_dlog_thalf <- t_half * dPhi_dt_half
  dphi_dlog_thalf <- t_half * dphi_dt_half

  # Derivative w.r.t. nu
  h_nu <- eps_rel * max(abs(nu), 1)
  d_plus  <- perturb_decompos(t_half, nu + h_nu, m)
  d_minus <- perturb_decompos(t_half, nu - h_nu, m)
  if (!is.null(d_plus) && !is.null(d_minus)) {
    e_plus  <- extract(d_plus, type)
    e_minus <- extract(d_minus, type)
    dPhi_dnu <- (e_plus$Phi - e_minus$Phi) / (2 * h_nu)
    dphi_dnu <- (e_plus$phi - e_minus$phi) / (2 * h_nu)
  } else if (!is.null(d_plus)) {
    e_plus <- extract(d_plus, type)
    dPhi_dnu <- (e_plus$Phi - base$Phi) / h_nu
    dphi_dnu <- (e_plus$phi - base$phi) / h_nu
  } else if (!is.null(d_minus)) {
    e_minus <- extract(d_minus, type)
    dPhi_dnu <- (base$Phi - e_minus$Phi) / h_nu
    dphi_dnu <- (base$phi - e_minus$phi) / h_nu
  } else {
    dPhi_dnu <- rep(0, n)
    dphi_dnu <- rep(0, n)
  }

  # Derivative w.r.t. m
  h_m <- eps_rel * max(abs(m), 1)
  d_plus  <- perturb_decompos(t_half, nu, m + h_m)
  d_minus <- perturb_decompos(t_half, nu, m - h_m)
  if (!is.null(d_plus) && !is.null(d_minus)) {
    e_plus  <- extract(d_plus, type)
    e_minus <- extract(d_minus, type)
    dPhi_dm <- (e_plus$Phi - e_minus$Phi) / (2 * h_m)
    dphi_dm <- (e_plus$phi - e_minus$phi) / (2 * h_m)
  } else if (!is.null(d_plus)) {
    e_plus <- extract(d_plus, type)
    dPhi_dm <- (e_plus$Phi - base$Phi) / h_m
    dphi_dm <- (e_plus$phi - base$phi) / h_m
  } else if (!is.null(d_minus)) {
    e_minus <- extract(d_minus, type)
    dPhi_dm <- (base$Phi - e_minus$Phi) / h_m
    dphi_dm <- (base$phi - e_minus$phi) / h_m
  } else {
    dPhi_dm <- rep(0, n)
    dphi_dm <- rep(0, n)
  }

  list(
    Phi             = base$Phi,
    phi             = base$phi,
    dPhi_dlog_thalf = dPhi_dlog_thalf,
    dPhi_dnu        = dPhi_dnu,
    dPhi_dm         = dPhi_dm,
    dphi_dlog_thalf = dphi_dlog_thalf,
    dphi_dnu        = dphi_dnu,
    dphi_dm         = dphi_dm
  )
}


#' Instantaneous hazard contribution from a single phase
#'
#' Computes \eqn{\phi_j(t) = d\Phi_j/dt} for one phase — the derivative of
#' the cumulative hazard contribution returned by [hzr_phase_cumhaz()].
#'
#' @inheritParams hzr_phase_cumhaz
#'
#' @return Numeric vector of instantaneous hazard contributions \eqn{\phi(t)},
#'   same length as `time`.
#'
#' @details
#' - `"cdf"`: \eqn{\phi(t) = g(t)} (density).
#' - `"hazard"`: \eqn{\phi(t) = h(t) = g(t)/(1-G(t))}.
#' - `"constant"`: \eqn{\phi(t) = 1}.
#'
#' @examples
#' t_grid <- seq(0.1, 10, by = 0.1)
#' phi_early <- hzr_phase_hazard(t_grid, t_half = 2, nu = 2, m = 0,
#'                                type = "cdf")
#' phi_late  <- hzr_phase_hazard(t_grid, t_half = 5, nu = 1, m = 0,
#'                                type = "hazard")
#'
#' @seealso [hzr_decompos()] for the underlying parametric family,
#'   [hzr_phase_cumhaz()] for the cumulative version.
#'
#' @export
hzr_phase_hazard <- function(time, t_half = 1, nu = 1, m = 0,
                              type = c("cdf", "hazard", "constant")) {
  type <- match.arg(type)

  if (type == "constant") {
    return(rep(1, length(time)))
  }

  d <- hzr_decompos(time, t_half = t_half, nu = nu, m = m)

  switch(type,
    cdf    = d$g,
    hazard = d$h
  )
}
