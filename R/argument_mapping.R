# argument_mapping.R — SAS HAZARD → R argument translation table
#
# PURPOSE
# -------
# During the SAS-to-R migration, analysts may need to know how legacy HAZARD
# procedure statements map to hazard() arguments.  hzr_argument_mapping()
# exposes the mapping as a data frame so it can be consulted programmatically
# or rendered in vignettes.
#
# The internal table .hzr_argument_mapping_table drives the exported function.
# Add rows as new arguments are implemented; update implementation_status from
# "planned" to "implemented" when the feature lands.

#' Legacy HAZARD to hvtiRhazard argument mapping
#'
#' Returns a formal mapping table that defines how legacy SAS HAZARD/C-style
#' inputs map to `hazard(...)` arguments in this package.
#'
#' @param include_planned Logical; if `FALSE`, only rows marked as implemented
#'   are returned.
#'
#' @return A data frame with one row per mapping rule.
#' @examples
#' hzr_argument_mapping()
#' hzr_argument_mapping(include_planned = FALSE)
#' @export
hzr_argument_mapping <- function(include_planned = TRUE) {
  map <- .hzr_argument_mapping_table
  if (!isTRUE(include_planned)) {
    map <- map[map$implementation_status == "implemented", , drop = FALSE]
  }
  rownames(map) <- NULL
  map
}

.hzr_argument_mapping_table <- data.frame(
  sas_statement = c(
    # --- Core hazard() arguments -----------------------------------------------
    "HAZARD", "HAZARD", "HAZARD", "HAZARD", "HAZARD", "HAZARD", "HAZARD",
    "TIME", "EVENT", "PARMS", "DIST",
    # --- Multiphase (phases / hzr_phase) arguments -----------------------------
    "HAZARD", "HAZARD",
    "G1", "G1", "G1", "G1",
    "G2",
    "G3", "G3", "G3", "G3"
  ),
  legacy_input = c(
    "TIME variable", "EVENT/censor variable", "X covariate block",
    "initial parameters", "baseline distribution", "control options",
    "additional legacy options", "t", "status", "theta0", "dist",
    # Multiphase
    "phases (3-phase structure)", "MU_1, MU_2, MU_3",
    "THALF / RHO (early)", "NU (early)", "M (early)", "DELTA (early)",
    "G2 constant phase",
    "TAU (late)", "GAMMA (late)", "ALPHA (late)", "ETA (late)"
  ),
  c_concept = c(
    "obs time array", "event indicator array", "design matrix",
    "parameter vector", "phase distribution selector", "optimizer/control struct",
    "misc legacy switches", "time vector", "event vector", "starting coef", "dist selector",
    # Multiphase
    "3-phase Early/Const/Late", "per-phase scale factors",
    "early half-life", "early time exponent", "early shape", "early time transform",
    "constant hazard rate phase",
    "late half-life", "late shape/exponent", "late shape", "late time exponent"
  ),
  r_parameter = c(
    "time", "status", "x", "theta", "dist", "control", "...",
    "time", "status", "theta", "dist",
    # Multiphase
    "phases (list of hzr_phase())", "mu (via exp(log_mu) in theta)",
    "hzr_phase(t_half=)", "hzr_phase(nu=)", "hzr_phase(m=)", "(absorbed by decompos)",
    "hzr_phase('constant')",
    "hzr_phase(t_half=)", "hzr_phase(nu=, m=)", "hzr_phase(m=)", "hzr_phase(nu=)"
  ),
  required = c(
    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
    # Multiphase
    FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE,
    FALSE,
    FALSE, FALSE, FALSE, FALSE
  ),
  expected_type = c(
    "numeric vector", "numeric/logical vector", "numeric matrix or data.frame",
    "numeric vector", "character scalar", "named list", "named arguments",
    "numeric vector", "numeric/logical vector", "numeric vector", "character scalar",
    # Multiphase
    "list of hzr_phase objects", "numeric (per-phase)",
    "positive scalar", "numeric scalar", "numeric scalar", "numeric scalar",
    "hzr_phase('constant')",
    "positive scalar", "numeric scalar", "numeric scalar", "numeric scalar"
  ),
  transform_rule = c(
    "pass through as numeric", "coerce to numeric 0/1", "data.frame -> data.matrix",
    "length must equal ncol(x) when x is present", "normalized lower-case label",
    "stored in spec$control", "stored in legacy_args for parity",
    "pass through", "coerce to numeric", "map PARMS/INITIAL to theta", "map DIST= to dist",
    # Multiphase
    "list(early=hzr_phase('cdf',...), constant=hzr_phase('constant'), late=hzr_phase('hazard',...))",
    "exp(alpha_j) in internal parameterization; estimated on log scale",
    "maps directly to hzr_phase(t_half=) starting value",
    "maps directly to hzr_phase(nu=) starting value",
    "maps directly to hzr_phase(m=) starting value",
    "time transform B(t) = (exp(delta*t)-1)/delta absorbed into decompos shape",
    "hzr_phase('constant') with no shape parameters",
    "maps to hzr_phase(t_half=) for late phase",
    "GAMMA collapses onto nu and m in decompos parameterization",
    "maps to hzr_phase(m=) for late phase",
    "maps to hzr_phase(nu=) for late phase"
  ),
  implementation_status = c(
    "implemented", "implemented", "implemented", "implemented", "implemented",
    "implemented", "implemented", "implemented", "implemented", "planned", "implemented",
    # Multiphase
    "implemented", "implemented",
    "implemented", "implemented", "implemented", "implemented",
    "implemented",
    "implemented", "implemented", "implemented", "implemented"
  ),
  notes = c(
    "Core observation time input.",
    "Event indicator currently retained as numeric in object$data$status.",
    "Future versions will support richer design encoding helpers.",
    "Used by predict.hazard as coefficient vector.",
    "Current default is 'weibull'; more options planned.",
    "Control list is stored and reserved for optimizer parity.",
    "Supports legacy-style pass-through options during migration.",
    "Canonical SAS migration uses TIME= mapping.",
    "Canonical SAS migration uses EVENT= mapping.",
    "SAS PARMS syntax parser not yet implemented.",
    "SAS DIST keyword maps directly to dist.",
    # Multiphase
    "Use dist='multiphase' with phases argument. N-phase generalization of legacy 3-phase model.",
    "Each phase has its own scale mu_j(x) = exp(alpha_j + x*beta_j). Starting value via hzr_phase().",
    "Half-life: time at which G(t_half) = 0.5. Same concept as SAS RHO/THALF.",
    "Time exponent controlling rate dynamics. Same parameter name as SAS early NU.",
    "Shape exponent controlling distributional form. Same parameter name as SAS early M.",
    "The C DELTA controlled B(t) = (exp(delta*t)-1)/delta. This transform is absorbed by decompos().",
    "Flat background rate. No shape parameters estimated. SAS G2 equivalent.",
    "Late-phase half-life. TAU maps to the t_half concept via the decompos parameterization.",
    "GAMMA jointly determines nu and m for the late phase through the decompos reparameterization.",
    "ALPHA maps to the shape exponent m in the late-phase decompos parameterization.",
    "ETA maps to the time exponent nu in the late-phase decompos parameterization.")
  ,
  stringsAsFactors = FALSE
)
