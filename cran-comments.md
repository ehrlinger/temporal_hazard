# CRAN submission comments -- TemporalHazard 1.0.1

## Summary

This is a resubmission of v1.0.0 addressing every item raised by
Benjamin Altmann's review of 2026-05-13. All three requested changes have
been applied:

1. **`\value` tags in Rd files.** Added `@return` documentation to all seven
   functions that were missing it: `hazard()`, `coef.hazard()`,
   `vcov.hazard()`, `print.hzr_calibrate()`, `print.hzr_deciles()`,
   `print.hzr_gof()`, and `print.hzr_kaplan()`. Each describes the class and
   structure of the return value (or documents invisible return for print
   methods).

2. **Writing to home filespace.** Internal fixture generators in
   `R/golden_fixtures.R` previously fell back to `"inst/fixtures"` (a path
   relative to `getwd()`) when the package was not installed. The fallback now
   uses `tempdir()` so no function writes to the user's home filespace by
   default.

3. **Setting a specific seed.** `set.seed(42)` has been removed as an
   unconditional call inside each generator. Generators now accept an optional
   `seed` argument; when supplied, the global `.Random.seed` is saved before
   and restored via `on.exit()` after the call, leaving the caller's RNG state
   unmodified.

## Test environments

* **Local:** R 4.6.0 on macOS (aarch64-apple-darwin23).
  `R CMD check --as-cran` returns 0 errors, 0 warnings, 0 notes.
* **GitHub Actions matrix:** ubuntu-latest (R-devel / R-release / R-oldrel-1),
  macos-latest (R-release), windows-latest (R-release) — all green.
* **Reverse-dependency check:** `tools::package_dependencies(reverse = TRUE)`
  returns 0.

## NOTE disposition

No notes in local `R CMD check --as-cran`.

## Background

* Package contains five clinical reference datasets (~120 KB compressed) used
  in vignettes and parity tests against the original C/SAS HAZARD program.
* Vignettes use Quarto (`VignetteBuilder: quarto`).
* A subset of tests (SAS parity, multi-start multiphase fits, OMC raw-file
  derivation) are guarded with `skip_on_cran()` to stay within CRAN runtime
  budgets; the full suite runs on GitHub Actions.
* DOI URLs in vignettes (doi.org) may return HTTP 403 to automated checkers
  but resolve correctly in browsers.
