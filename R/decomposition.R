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
