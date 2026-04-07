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
    "HAZARD", "HAZARD", "HAZARD", "HAZARD", "HAZARD", "HAZARD", "HAZARD",
    "TIME", "EVENT", "PARMS", "DIST"
  ),
  legacy_input = c(
    "TIME variable", "EVENT/censor variable", "X covariate block",
    "initial parameters", "baseline distribution", "control options",
    "additional legacy options", "t", "status", "theta0", "dist"
  ),
  c_concept = c(
    "obs time array", "event indicator array", "design matrix",
    "parameter vector", "phase distribution selector", "optimizer/control struct",
    "misc legacy switches", "time vector", "event vector", "starting coef", "dist selector"
  ),
  r_parameter = c(
    "time", "status", "x", "theta", "dist", "control", "...",
    "time", "status", "theta", "dist"
  ),
  required = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  expected_type = c(
    "numeric vector", "numeric/logical vector", "numeric matrix or data.frame",
    "numeric vector", "character scalar", "named list", "named arguments",
    "numeric vector", "numeric/logical vector", "numeric vector", "character scalar"
  ),
  transform_rule = c(
    "pass through as numeric", "coerce to numeric 0/1", "data.frame -> data.matrix",
    "length must equal ncol(x) when x is present", "normalized lower-case label",
    "stored in spec$control", "stored in legacy_args for parity",
    "pass through", "coerce to numeric", "map PARMS/INITIAL to theta", "map DIST= to dist"
  ),
  implementation_status = c(
    "implemented", "implemented", "implemented", "implemented", "implemented",
    "implemented", "implemented", "implemented", "implemented", "planned", "implemented"
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
    "SAS DIST keyword maps directly to dist.")
  ,
  stringsAsFactors = FALSE
)
