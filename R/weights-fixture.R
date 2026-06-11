# weights-fixture.R -- Load and validate SAS fractional-weight parity fixtures.
#
# Fixtures live under `inst/fixtures/weighted-<scenario>.rds`.  Like the
# stepwise fixtures they are optional: the parity test skips when the file
# is absent, so the package installs and CI passes without any SAS output.
#
# The workflow for capturing one is documented in
# `inst/extdata/weights-fixtures/README.md`.  .hzr_build_weighted_fixture()
# packages the SAS coefficient table + meta into the expected list shape;
# .hzr_load_weighted_fixture() is the lookup called from the parity test.
#
# `.hzr_parse_keyvals()` is shared with stepwise-fixture.R.

#' Required top-level fields of a weighted-fit parity fixture
#'
#' @keywords internal
#' @noRd
.hzr_weighted_fixture_schema <- function() {
  list(
    meta = c(
      "dataset", "dist", "weight_var", "weight_kind",
      "sas_version", "captured_on", "proc_hazard"
    ),
    final = c("coef", "logLik")
  )
}


#' Validate a weighted-fit parity fixture against the schema
#'
#' @param fix A list as produced by `.hzr_build_weighted_fixture()` or
#'   deserialised from `inst/fixtures/weighted-*.rds`.
#' @return The fixture, unchanged, if valid.  Throws otherwise.
#'
#' @keywords internal
#' @noRd
.hzr_validate_weighted_fixture <- function(fix) {
  schema <- .hzr_weighted_fixture_schema()
  problems <- character()

  for (group in names(schema)) {
    if (!group %in% names(fix)) {
      problems <- c(problems, paste0("missing top-level: ", group))
      next
    }
    missing_fields <- setdiff(schema[[group]], names(fix[[group]]))
    if (length(missing_fields) > 0L) {
      problems <- c(problems,
                    sprintf("missing %s: %s", group,
                            paste(missing_fields, collapse = ", ")))
    }
  }

  if (!is.null(fix$final) && !is.null(fix$final$coef)) {
    needed <- c("variable", "phase", "estimate")
    miss   <- setdiff(needed, names(fix$final$coef))
    if (length(miss) > 0L) {
      problems <- c(problems,
                    sprintf("final$coef missing columns: %s",
                            paste(miss, collapse = ", ")))
    }
  }

  if (length(problems) > 0L) {
    stop("Invalid weighted parity fixture:\n  - ",
         paste(problems, collapse = "\n  - "),
         call. = FALSE)
  }
  invisible(fix)
}


#' Load a weighted-fit parity fixture by name
#'
#' Returns NULL (not an error) when the fixture file is absent.
#'
#' @param name Short name like `"avc-fractional"`.
#' @return The fixture list or NULL.
#'
#' @keywords internal
#' @noRd
.hzr_load_weighted_fixture <- function(name) {
  file_name <- paste0("weighted-", name, ".rds")
  path <- system.file("fixtures", file_name, package = "TemporalHazard")
  if (!nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  fix <- readRDS(path)
  .hzr_validate_weighted_fixture(fix)
}


#' Build a weighted-fit parity fixture from SAS capture files
#'
#' Reads the coefficient CSV and `key=value` meta file emitted by
#' `inst/extdata/weights-fixtures/parse-weighted-lst.R` and assembles the
#' fixture list.  Writes the result to `out_path` if supplied.
#'
#' @param coef_csv Path to the parsed coefficient table CSV (columns
#'   `variable`, `phase`, `estimate`, and optionally `std_error`,
#'   `z_stat`, `p_value`).
#' @param meta_txt Path to the `key=value` meta file (must include a
#'   `logLik=` line, added by the parser).
#' @param out_path Optional output `.rds` path.
#'
#' @return The fixture list.
#'
#' @keywords internal
#' @noRd
.hzr_build_weighted_fixture <- function(coef_csv, meta_txt, out_path = NULL) {
  if (!file.exists(coef_csv)) {
    stop("coef_csv not found: ", coef_csv, call. = FALSE)
  }
  if (!file.exists(meta_txt)) {
    stop("meta_txt not found: ", meta_txt, call. = FALSE)
  }

  coef_tbl <- utils::read.csv(coef_csv, stringsAsFactors = FALSE)
  if ("phase" %in% names(coef_tbl)) {
    coef_tbl$phase[coef_tbl$phase == "."] <- NA_character_
  }

  meta_kv <- .hzr_parse_keyvals(readLines(meta_txt))

  fix <- list(
    meta = list(
      dataset     = meta_kv[["dataset"]]     %||% NA_character_,
      dist        = meta_kv[["dist"]]        %||% "multiphase",
      model       = meta_kv[["model"]]       %||% NA_character_,
      covariates  = meta_kv[["covariates"]]  %||% NA_character_,
      weight_var  = meta_kv[["weight_var"]]  %||% NA_character_,
      weight_kind = meta_kv[["weight_kind"]] %||% NA_character_,
      conserve    = as.integer(meta_kv[["conserve"]] %||% NA_integer_),
      sas_version = meta_kv[["sas_version"]] %||% NA_character_,
      captured_on = meta_kv[["captured_on"]] %||% NA_character_,
      proc_hazard = meta_kv[["proc_hazard"]] %||% NA_character_,
      notes       = meta_kv[["notes"]]
    ),
    final = list(
      coef   = coef_tbl,
      logLik = as.numeric(meta_kv[["logLik"]] %||% NA_real_)
    )
  )

  .hzr_validate_weighted_fixture(fix)

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(fix, out_path)
  }
  fix
}
