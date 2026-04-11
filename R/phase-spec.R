# phase-spec.R — Phase specification for multiphase hazard models
#
# PURPOSE
# -------
# Defines the hzr_phase() constructor that specifies one phase in an
# N-phase additive cumulative hazard model:
#
#   H(t | x) = sum_j  mu_j(x) * Phi_j(t; t_half_j, nu_j, m_j)
#
# Each phase is one term in the sum.  The constructor stores:
#   - type:    "cdf" (early risk), "hazard" (late risk), or "constant"
#   - t_half, nu, m: starting values for the decompos shape parameters
#   - formula: optional one-sided formula for phase-specific covariates
#
# The helpers extract metadata needed during likelihood construction:
#   - .hzr_phase_n_shape():     number of shape parameters (3 or 0)
#   - .hzr_phase_theta_names(): named labels for the parameter sub-vector
#
# SAS/C BRIDGE
# ------------
# SAS users specify Early/Constant/Late phases via DIST keywords.  The
# hzr_phase() constructor maps directly:
#   SAS Early (G1) -> hzr_phase("cdf",    t_half, nu, m)
#   SAS Const (G2) -> hzr_phase("constant")
#   SAS Late  (G3) -> hzr_phase("hazard", t_half, nu, m)

# ============================================================================
# Constructor
# ============================================================================

#' Specify a single hazard phase
#'
#' Creates an `hzr_phase` object describing one term in a multiphase additive
#' cumulative hazard model.  Pass a list of these to the `phases` argument of
#' [hazard()] when `dist = "multiphase"`.
#'
#' @section Phase types:
#'
#' \describe{
#'   \item{`"cdf"`}{Early risk that resolves over time.
#'     \eqn{\Phi(t) = G(t)}, bounded \eqn{[0, 1]}.
#'     SAS equivalent: Early / G1 phase.}
#'   \item{`"hazard"`}{Late or aging risk that accumulates.
#'     \eqn{\Phi(t) = -\log(1 - G(t))}, monotone increasing.
#'     SAS equivalent: Late / G3 phase.}
#'   \item{`"constant"`}{Flat background hazard rate.
#'     \eqn{\Phi(t) = t}. No shape parameters are estimated.
#'     SAS equivalent: Constant / G2 phase.}
#' }
#'
#' @param type Character; phase type: `"cdf"`, `"hazard"`, or `"constant"`.
#' @param t_half Positive scalar; initial half-life (time at which
#'   \eqn{G(t_{1/2}) = 0.5}).  Ignored for `"constant"` phases.
#'   SAS early: `THALF`/`RHO`.  SAS late: relates to `TAU`.
#' @param nu Numeric scalar; initial time exponent.  Ignored for `"constant"`.
#'   SAS early: `NU`.  SAS late: relates to `GAMMA`/`ETA`.
#' @param m Numeric scalar; initial shape exponent.  Ignored for `"constant"`.
#'   SAS early: `M`.  SAS late: relates to `GAMMA`/`ALPHA`.
#' @param formula Optional one-sided formula (e.g. `~ age + nyha`) for
#'   phase-specific covariates.  When `NULL` (default), the phase inherits
#'   the global formula from [hazard()].
#'
#' @return An S3 object of class `"hzr_phase"` with elements:
#' \describe{
#'   \item{type}{Phase type string.}
#'   \item{t_half}{Initial half-life.}
#'   \item{nu}{Initial time exponent.}
#'   \item{m}{Initial shape exponent.}
#'   \item{formula}{Phase-specific formula or `NULL`.}
#' }
#'
#' @examples
#' # Classic 3-phase HAZARD pattern
#' early <- hzr_phase("cdf",      t_half = 0.5, nu = 2, m = 0)
#' const <- hzr_phase("constant")
#' late  <- hzr_phase("hazard",   t_half = 5,   nu = 1, m = 0)
#'
#' # Phase with specific covariates
#' early_cov <- hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
#'                         formula = ~ age + shock)
#'
#' # Use in hazard():
#' # hazard(Surv(time, status) ~ age, data = dat,
#' #        dist = "multiphase",
#' #        phases = list(early = early, constant = const, late = late))
#'
#' @seealso [hzr_decompos()] for the underlying parametric family,
#'   [hzr_phase_cumhaz()] and [hzr_phase_hazard()] for computing
#'   \eqn{\Phi(t)} and \eqn{\phi(t)} from these specifications.
#'
#' @export
hzr_phase <- function(type = c("cdf", "hazard", "constant"),
                      t_half = 1, nu = 1, m = 0,
                      formula = NULL) {

  type <- match.arg(type)

  # --- Validate shape parameters (non-constant phases only) -----------------
  if (type != "constant") {
    stopifnot(
      "t_half must be a positive scalar" =
        is.numeric(t_half) && length(t_half) == 1L &&
        is.finite(t_half) && t_half > 0,
      "nu must be a finite numeric scalar" =
        is.numeric(nu) && length(nu) == 1L && is.finite(nu),
      "m must be a finite numeric scalar" =
        is.numeric(m) && length(m) == 1L && is.finite(m)
    )

    if (m < 0 && nu < 0) {
      stop("Decomposition undefined when both m < 0 and nu < 0 (m = ",
           m, ", nu = ", nu, ").  See ?hzr_decompos for valid cases.",
           call. = FALSE)
    }
  }

  # --- Validate formula if supplied -----------------------------------------
  if (!is.null(formula)) {
    if (!inherits(formula, "formula")) {
      stop("formula must be a one-sided formula (e.g. ~ age + nyha) or NULL.",
           call. = FALSE)
    }
    # Must be one-sided (no LHS)
    if (length(formula) == 3L) {
      stop("Phase formula must be one-sided (e.g. ~ age + nyha), ",
           "not two-sided (response ~ predictors).", call. = FALSE)
    }
  }

  structure(
    list(
      type    = type,
      t_half  = if (type == "constant") NA_real_ else t_half,
      nu      = if (type == "constant") NA_real_ else nu,
      m       = if (type == "constant") NA_real_ else m,
      formula = formula
    ),
    class = "hzr_phase"
  )
}


# ============================================================================
# S3 methods
# ============================================================================

#' @export
print.hzr_phase <- function(x, ...) {
  label <- switch(x$type,
    cdf      = "cdf (early risk)",
    hazard   = "hazard (late risk)",
    constant = "constant (flat rate)"
  )
  cat("<hzr_phase>", label, "\n")

  if (x$type != "constant") {
    cat("  t_half =", format(x$t_half, digits = 4),
        " nu =", format(x$nu, digits = 4),
        " m =", format(x$m, digits = 4), "\n")
  }

  if (!is.null(x$formula)) {
    cat("  covariates:", deparse(x$formula), "\n")
  }

  invisible(x)
}


# ============================================================================
# Type-checking helper
# ============================================================================

#' Test if an object is an hzr_phase
#'
#' @param x Object to test.
#' @return Logical scalar.
#'
#' @examples
#' is_hzr_phase(hzr_phase("cdf"))
#' is_hzr_phase("not a phase")
#'
#' @export
is_hzr_phase <- function(x) {
  inherits(x, "hzr_phase")
}


# ============================================================================
# Internal helpers for likelihood construction
# ============================================================================

#' Number of shape parameters for a phase
#'
#' Returns 3 (t_half, nu, m) for `"cdf"` and `"hazard"` phases, 0 for
#' `"constant"`.
#'
#' @param phase An `hzr_phase` object.
#' @return Integer: 3 or 0.
#' @keywords internal
.hzr_phase_n_shape <- function(phase) {
  stopifnot(is_hzr_phase(phase))
  if (phase$type == "constant") 0L else 3L
}


#' Generate parameter names for a phase's sub-vector
#'
#' Builds the named labels for the parameter block that belongs to a single
#' phase in the full theta vector.  The block layout is:
#'
#' \itemize{
#'   \item For `"cdf"`/`"hazard"`: `[log_mu, log_t_half, nu, m, beta_1, ..., beta_p]`
#'   \item For `"constant"`:       `[log_mu, beta_1, ..., beta_p]`
#' }
#'
#' @param phase An `hzr_phase` object.
#' @param phase_name Character label for the phase (e.g. `"early"`, `"phase_1"`).
#' @param covariate_names Character vector of covariate column names that
#'   this phase uses.  Can be length 0 if no covariates.
#' @return Character vector of parameter names.
#' @keywords internal
.hzr_phase_theta_names <- function(phase, phase_name, covariate_names = character(0)) {
  stopifnot(is_hzr_phase(phase))

  # Intercept (always present)
  names_out <- paste0(phase_name, ".log_mu")

  # Shape parameters (cdf/hazard only)
  if (phase$type != "constant") {
    names_out <- c(names_out,
      paste0(phase_name, ".log_t_half"),
      paste0(phase_name, ".nu"),
      paste0(phase_name, ".m")
    )
  }

  # Covariate coefficients
  if (length(covariate_names) > 0L) {
    names_out <- c(names_out,
      paste0(phase_name, ".", covariate_names)
    )
  }

  names_out
}


#' Total number of parameters for a phase
#'
#' Returns the count of free parameters in the theta sub-vector for one phase:
#' 1 (log_mu) + n_shape + n_covariates.
#'
#' @param phase An `hzr_phase` object.
#' @param n_covariates Integer; number of covariate columns this phase uses.
#' @return Integer.
#' @keywords internal
.hzr_phase_n_params <- function(phase, n_covariates = 0L) {
  stopifnot(is_hzr_phase(phase))
  1L + .hzr_phase_n_shape(phase) + as.integer(n_covariates)
}


#' Extract starting values from a phase specification
#'
#' Returns initial theta sub-vector on the estimation (internal) scale:
#' log(mu), log(t_half), nu, m, followed by zeros for covariate coefficients.
#'
#' @param phase An `hzr_phase` object.
#' @param n_covariates Integer; number of covariate columns.
#' @param mu_start Numeric scalar; initial scale parameter (default 0.1).
#' @return Named numeric vector of starting values.
#' @keywords internal
.hzr_phase_start <- function(phase, n_covariates = 0L, mu_start = 0.1) {
  stopifnot(is_hzr_phase(phase))
  stopifnot(
    "mu_start must be a positive scalar" =
      is.numeric(mu_start) && length(mu_start) == 1L && mu_start > 0
  )

  vals <- log(mu_start)  # log_mu

  if (phase$type != "constant") {
    vals <- c(vals, log(phase$t_half), phase$nu, phase$m)
  }

  # Covariate betas initialized to zero
  if (n_covariates > 0L) {
    vals <- c(vals, rep(0, n_covariates))
  }

  vals
}


#' Validate a list of phase specifications
#'
#' Checks that `phases` is a non-empty named list of `hzr_phase` objects.
#' Auto-names unnamed phases as `phase_1`, `phase_2`, etc.
#'
#' @param phases A list of `hzr_phase` objects.
#' @return The validated (and possibly auto-named) list, invisibly.
#' @keywords internal
.hzr_validate_phases <- function(phases) {
  # Catch bare hzr_phase passed instead of list(hzr_phase(...))
  if (is_hzr_phase(phases)) {
    stop("phases must be a non-empty list of hzr_phase objects, ",
         "not a single hzr_phase.  Wrap it: list(hzr_phase(...)).", call. = FALSE)
  }
  if (!is.list(phases) || length(phases) == 0L) {
    stop("phases must be a non-empty list of hzr_phase objects.", call. = FALSE)
  }

  # Check each element

  for (i in seq_along(phases)) {
    if (!is_hzr_phase(phases[[i]])) {
      stop("phases[[", i, "]] is not an hzr_phase object. ",
           "Use hzr_phase() to create phase specifications.", call. = FALSE)
    }
  }

  # Auto-name unnamed phases
  nms <- names(phases)
  if (is.null(nms)) {
    nms <- character(length(phases))
  }
  empty <- nms == "" | is.na(nms)
  if (any(empty)) {
    nms[empty] <- paste0("phase_", which(empty))
  }

  # Check for duplicate names
  if (anyDuplicated(nms)) {
    stop("Phase names must be unique. Duplicates found: ",
         paste(nms[duplicated(nms)], collapse = ", "), call. = FALSE)
  }

  names(phases) <- nms
  invisible(phases)
}
