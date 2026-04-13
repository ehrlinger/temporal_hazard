## ── Convert inst/extdata CSVs to lazy-loaded data/ objects ──────────────────
##
## Run this script from the package root:
##   source("data-raw/make_data.R")
##
## It reads each CSV from inst/extdata/, creates the corresponding data.frame,
## and saves it as a compressed .rda file in data/.
##
## All CSVs originate from SAS exports where "." denotes missing values.

dir.create("data", showWarnings = FALSE)

.read_sas_csv <- function(path, ...) {
  read.csv(path, stringsAsFactors = FALSE, na.strings = c("NA", "."), ...)
}

# ── AVC: Aortic Valve Replacement (Cleveland Clinic, 1977–1993) ─────────────
avc <- .read_sas_csv("inst/extdata/avc.csv")
stopifnot(nrow(avc) == 310)
usethis::use_data(avc, overwrite = TRUE, compress = "xz")

# ── CABGKUL: Primary Isolated CABG (KU Leuven, 1971–2008) ──────────────────
cabgkul <- .read_sas_csv("inst/extdata/cabgkul.csv")
stopifnot(nrow(cabgkul) == 5880, sum(cabgkul$dead) == 545)
usethis::use_data(cabgkul, overwrite = TRUE, compress = "xz")

# ── OMC: Open Mitral Commissurotomy ─────────────────────────────────────────
omc <- .read_sas_csv("inst/extdata/omc.csv")
stopifnot(nrow(omc) == 339)
usethis::use_data(omc, overwrite = TRUE, compress = "xz")

# ── TGA: Transposition of the Great Arteries (BCH/CHOP) ────────────────────
tga <- .read_sas_csv("inst/extdata/tga.csv")
stopifnot(nrow(tga) == 470)
usethis::use_data(tga, overwrite = TRUE, compress = "xz")

# ── Valves: Heart Valve Replacement ─────────────────────────────────────────
valves <- .read_sas_csv("inst/extdata/valves.csv")
stopifnot(nrow(valves) == 1533)
usethis::use_data(valves, overwrite = TRUE, compress = "xz")

message("All 5 datasets saved to data/")
