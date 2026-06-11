# make-avc-weighted-csv.R -- regenerate payload/avc-weighted.csv.
#
# Writes the AVC analysis data with a deterministic, NON-INTEGER weight
# column `ipw` appended.  The weight is a closed-form function of `age`
# only -- no RNG, no cohort-derived statistics -- so the column is exactly
# reproducible and SAS (PROC HAZARD `WEIGHT ipw;`) and R (`hazard(weights =
# ipw)`) read identical values.  The precise formula is illustrative (an
# inverse-probability-style weight in (0.5, 1.5)); what matters for the
# parity test is only that both engines weight the same fractional column.
#
# Run from the project root after devtools::load_all():
#   source("inst/extdata/weights-fixtures/make-avc-weighted-csv.R")

suppressMessages(devtools::load_all(quiet = TRUE))
data(avc, package = "TemporalHazard")

# Deterministic fractional weight: logistic in age, centred at 60, scale 10.
# Range (0.5, 1.5); rounded to 6 dp so the CSV is byte-stable.
avc$ipw <- round(0.5 + 1 / (1 + exp(-(avc$age - 60) / 10)), 6)

out <- "inst/extdata/weights-fixtures/payload/avc-weighted.csv"
utils::write.csv(avc, out, row.names = FALSE, na = ".")
message("Wrote ", out, " (", nrow(avc), " rows; ipw range ",
        paste(round(range(avc$ipw, na.rm = TRUE), 4), collapse = "-"), ")")
