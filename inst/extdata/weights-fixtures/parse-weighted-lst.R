# parse-weighted-lst.R -- Parse a HAZARD .lst listing from the weighted fit
# into the files consumed by .hzr_build_weighted_fixture().
#
# Usage (from the project root, after the SAS run via payload/run.sh):
#
#   source("inst/extdata/weights-fixtures/parse-weighted-lst.R")
#
# Inputs:
#   inst/extdata/weights-fixtures/payload/capture/avc-weighted.lst
#   inst/extdata/weights-fixtures/payload/capture/weighted_meta.txt
#
# Outputs (written in-place to payload/capture/):
#   weighted_coef.csv  -- final coefficient table
#   weighted_meta.txt  -- updated with logLik= entry
#
# Then builds:
#   inst/fixtures/weighted-avc-fractional.rds
#
# The coefficient-table and log-likelihood parsing mirror
# stepwise-fixtures/parse-hazard-lst.R (same PROC HAZARD listing format).

lst_file    <- "inst/extdata/weights-fixtures/payload/capture/avc-weighted.lst"
meta_file   <- "inst/extdata/weights-fixtures/payload/capture/weighted_meta.txt"
out_coef    <- "inst/extdata/weights-fixtures/payload/capture/weighted_coef.csv"
fixture_out <- "inst/fixtures/weighted-avc-fractional.rds"

stopifnot(file.exists(lst_file), file.exists(meta_file))

lines <- readLines(lst_file, warn = FALSE)

# ---- Final coefficient table ------------------------------------------------
# Lives in the "Parameter Estimate Summary" block.  Fixed-width columns:
#
#   Phase      Parameter    Estimate       Std error             Z    Prob>|Z|
#   ----------------------...
#   Early:     E0             -2.20788      0.4338312         -5.089    <.0001
#              AGE          -0.00544498    0.002572203         -2.117    0.0343
#              ------...
#   Constant:  C0             -10.6589       1.726197         -6.175    <.0001
#              AGE           ...
#   ------...

param_summary_idx <- grep("Parameter Estimate Summary", lines, fixed = TRUE)[1]
stopifnot(!is.na(param_summary_idx))

col_header_idx <- grep(
  "Phase\\s+Parameter\\s+Estimate\\s+Std error",
  lines[param_summary_idx:length(lines)]
)[1] + param_summary_idx - 1L

data_start <- col_header_idx + 2L   # +1 dashes, +1 first data line

stop_idx <- grep("Estimates for Model Parameters", lines, fixed = TRUE)
stop_idx <- stop_idx[stop_idx > data_start][1L]

final_rows    <- list()
current_phase <- NA_character_

for (i in seq(data_start, stop_idx - 1L)) {
  line    <- lines[i]
  trimmed <- trimws(line)

  if (nchar(trimmed) == 0L) next
  if (grepl("^-{5,}", trimmed))  next

  phase_lbl <- regmatches(
    line, regexpr("(Early|Constant):", line, perl = TRUE, ignore.case = TRUE)
  )
  if (length(phase_lbl) > 0L && nchar(phase_lbl) > 0L) {
    current_phase <- tolower(sub(":", "", phase_lbl))
    line <- sub(".*(?:Early|Constant):\\s+", "", line, ignore.case = TRUE,
                perl = TRUE)
  }

  m <- regexpr(
    paste0("^\\s*(\\S+)\\s+([+-]?[0-9.Ee+-]+)\\s+([0-9.Ee+-]+)\\s+",
           "([+-]?[0-9.Ee+-]+)\\s+([<]?[.0-9]+)"),
    line, perl = TRUE
  )
  if (m < 0L) next

  tok <- strsplit(trimws(regmatches(line, m)), "\\s+")[[1]]
  if (length(tok) < 5L) next

  p_raw <- tok[5L]
  p_val <- if (startsWith(p_raw, "<")) {
    as.numeric(sub("^<\\.?", "0.", p_raw))
  } else {
    as.numeric(p_raw)
  }

  final_rows[[length(final_rows) + 1L]] <- data.frame(
    variable  = tok[1L],
    phase     = current_phase,
    estimate  = as.numeric(tok[2L]),
    std_error = as.numeric(tok[3L]),
    z_stat    = as.numeric(tok[4L]),
    p_value   = p_val,
    stringsAsFactors = FALSE
  )
}

final_df <- do.call(rbind, final_rows)
final_df$variable <- tolower(final_df$variable)
utils::write.csv(final_df, out_coef, row.names = FALSE)
message("Wrote ", out_coef, " (", nrow(final_df), " parameters)")

# ---- Augment meta.txt with logLik -------------------------------------------
ll_val <- as.numeric(trimws(sub(
  ".*Log likelihood\\s*=\\s*", "",
  tail(grep("Log likelihood\\s*=", lines, value = TRUE), 1L)
)))
message("Final log-likelihood: ", ll_val)

meta_lines <- readLines(meta_file)
meta_lines <- meta_lines[!grepl("^logLik\\s*=", meta_lines)]
meta_lines <- c(meta_lines, paste0("logLik=", ll_val))
writeLines(meta_lines, meta_file)
message("Updated ", meta_file)

# ---- Build the .rds fixture -------------------------------------------------
devtools::load_all(quiet = TRUE)

fix <- TemporalHazard:::.hzr_build_weighted_fixture(
  coef_csv = out_coef,
  meta_txt = meta_file,
  out_path = fixture_out
)
message("Built fixture: ", fixture_out)
message("\nFinal coefficients:")
print(fix$final$coef)
message("\nlogLik: ", fix$final$logLik)
