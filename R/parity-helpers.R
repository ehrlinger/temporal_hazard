# parity-helpers.R — Helpers for cross-validating R output against the legacy C binary
#
# STATUS: Stub / scaffold (M1 infrastructure, awaiting C binary integration)
#
# PURPOSE
# -------
# These functions wrap the compiled HAZARD C binary so the test suite can invoke
# it programmatically, capture its output, and compare parameter estimates with
# those produced by the R implementation.  This cross-validation is the primary
# correctness signal during the SAS→R migration.
#
# CURRENT STATE
# -------------
# .hzr_get_hazard_binary() — locates the binary at inst/bin/hazard
# .hzr_run_hazard_binary() — scaffolded; writes input CSV and builds the command
#                             string, but the actual system() invocation and output
#                             parser are not yet finalised (see TODO comments).
# .hzr_parse_hazard_output() — parser stub; returns list(NA) until implemented.
#
# ENABLING C PARITY TESTS
# -----------------------
# 1. Place the compiled binary at inst/bin/hazard (chmod +x).
# 2. Complete .hzr_run_hazard_binary(): actual system() call + output parser.
# 3. Update test-parity-core.R to call .hzr_run_hazard_binary() instead of
#    loading .rds fixtures where true C-binary comparison is desired.
#
# Until the binary is wired up, parity tests use synthetic R-generated fixtures
# produced by golden_fixtures.R.

#' Parity Testing Helpers
#'
#' Functions to generate golden fixtures from the C reference \code{hazard} binary
#' and validate R implementation against them.
#'
#' @keywords internal

#' Get path to compiled hazard binary
#'
#' @noRd
.hzr_get_hazard_binary <- function() {
  # Path to compiled hazard executable in package inst/bin
  bin_path <- system.file("bin", "hazard", package = "TemporalHazard")
  if (!file.exists(bin_path)) {
    stop(
      "Hazard C binary not found at: ", bin_path, "\n",
      "Run: cp /path/to/hazard/src/hazard/hazard inst/bin/",
      call. = FALSE
    )
  }
  bin_path
}

#' Run hazard C binary on input data
#'
#' Invokes the compiled \code{hazard} executable with prepared input and parses output.
#'
#' @param time Numeric vector of follow-up times
#' @param status Numeric vector of event indicators (0/1)
#' @param x Optional design matrix (univariable if NULL; multivar if provided)
#' @param theta Starting parameter vector (names must match C code expectations)
#' @param dist Distribution name ("weibull", etc.)
#' @param control Optional control list (nocov, nocor, condition, etc.)
#' @param data_name Character name for temporary data file prefix
#'
#' @return List with:
#'   - \code{call_string}: Full command-line invocation for reproducibility
#'   - \code{stdout}: Raw standard output from hazard binary
#'   - \code{parsed}: Parsed results (tibble with estimates, SE, z, p-values)
#'   - \code{info}: Metadata (df, condition, convergence flags)
#'
#' @noRd
.hzr_run_hazard_binary <- function(
    time,
    status,
    x = NULL,
    theta = NULL,
    dist = "weibull",
    control = list(),
    data_name = "golden_test") {
  
  bin_path <- .hzr_get_hazard_binary()
  temp_dir <- tempdir()
  
  # Prepare data file in XPORT format (SAS transport, as hazard expects)
  data_file <- file.path(temp_dir, paste0(data_name, ".xpt"))
  
  # For now, write simple CSV and note that real implementation would use haven::write_xpt
  # or similar to create SAS XPORT format
  df <- data.frame(
    TIME = time,
    STATUS = as.integer(status)
  )
  
  if (!is.null(x)) {
    if (is.vector(x)) {
      df[["X1"]] <- x
    } else {
      for (i in seq_len(ncol(x))) {
        df[[paste0("X", i)]] <- x[, i]
      }
    }
  }
  
  # Write CSV (hazard binary expects specific format; this is a placeholder)
  write.csv(df, data_file, row.names = FALSE)
  
  # Build command line arguments
  # Format: hazard [-options] -indata=file -outdata=file ...
  cmd_args <- c(
    "-indata=" %+% data_file,
    "-outdata=" %+% file.path(temp_dir, paste0(data_name, "_out.txt"))
  )
  
  if (control$nocov %||% FALSE) cmd_args <- c(cmd_args, "-nocov")
  if (control$nocor %||% FALSE) cmd_args <- c(cmd_args, "-nocor")
  
  # Full command for logging
  call_string <- paste(bin_path, paste(cmd_args, collapse = " "))
  
  # Execute binary
  result <- tryCatch(
    system2(bin_path, args = cmd_args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      stop("Failed to run hazard binary: ", e$message, call. = FALSE)
    }
  )
  
  out_file <- file.path(temp_dir, paste0(data_name, "_out.txt"))
  out_text <- if (file.exists(out_file)) readLines(out_file) else ""
  
  # Parse output (to be implemented based on actual hazard output format)
  parsed <- .hzr_parse_hazard_output(out_text, theta = theta)
  
  list(
    call_string = call_string,
    stdout = paste(result, collapse = "\n"),
    parsed = parsed,
    info = list(
      dist = dist,
      n = length(time),
      events = sum(status),
      binary_path = bin_path
    )
  )
}

#' Parse hazard binary output (placeholder)
#'
#' TODO: Implement once the exact column layout of the C binary's output file
#' is documented.  Expected columns: parameter name, estimate, SE, z-statistic,
#' p-value.  See inst/bin/hazard --help for format details.
#'
#' @noRd
.hzr_parse_hazard_output <- function(output_text, theta = NULL) {
  # This function will be filled in once we understand
  # the exact output format of the compiled binary.
  # For now, return empty tibble with required columns.
  
  structure(
    tibble::tibble(
      parameter = character(),
      estimate = numeric(),
      std_err = numeric(),
      z_stat = numeric(),
      p_value = numeric()
    ),
    output_text = paste(output_text, collapse = "\n")
  )
}

#' Generate golden fixture for univariable Weibull model
#'
#' Runs the C hazard binary on a dataset and stores results as test fixture.
#'
#' @param fixture_name Character identifying the fixture (e.g., "avc_death_hz")
#' @param time Numeric vector of follow-up times
#' @param status Numeric vector of event indicators
#' @param x Optional design matrix for multivariable models
#' @param theta Optional starting parameters
#' @param control Optional control list
#' @param output_dir Directory to store fixture (default: \code{inst/fixtures/})
#'
#' @return Invisibly returns the parsed output of the reference binary.
#'
#' @noRd
.hzr_generate_golden_fixture <- function(
    fixture_name,
    time,
    status,
    x = NULL,
    theta = NULL,
    control = list(),
    output_dir = NULL) {
  
  if (is.null(output_dir)) {
    output_dir <- system.file("fixtures", package = "TemporalHazard")
    if (!nzchar(output_dir)) {
      output_dir <- "inst/fixtures"
    }
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Run binary
  result <- .hzr_run_hazard_binary(
    time = time,
    status = status,
    x = x,
    theta = theta,
    control = control,
    data_name = fixture_name
  )
  
  # Save fixture as RDS for reproducible testing
  fixture_file <- file.path(output_dir, paste0(fixture_name, ".rds"))
  saveRDS(result, fixture_file)
  
  message("Generated golden fixture: ", fixture_file)
  
  invisible(result$parsed)
}

#' Utility: string concatenation pipe
#'
#' @noRd
'%+%' <- function(x, y) paste0(x, y)

#' Utility: logical OR with NA handling
#'
#' @noRd
'%||%' <- function(x, y) if (is.null(x) || is.na(x)) y else x
