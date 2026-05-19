# CRAN submission comments -- TemporalHazard 1.0.2

## Summary

This is a resubmission addressing Benjamin Altmann's follow-up review. Two of
the three items from the prior review were resolved in 1.0.1; the
home-filespace item was only partially fixed and is now fully resolved in
1.0.2. Each point is itemised below.

1. **`\value` tags in Rd files.** Resolved (already in place since 1.0.1).
   `@return` documentation is present for all seven flagged functions —
   `hazard()`, `coef.hazard()`, `vcov.hazard()`, `print.hzr_calibrate()`,
   `print.hzr_deciles()`, `print.hzr_gof()`, `print.hzr_kaplan()` — each
   describing the class and structure of the return value. The four `print.*`
   methods document that they return their argument `x` invisibly (called for
   the printing side effect) and describe the columns/attributes of that
   object.

2. **Writing files to the home filespace.** Now fully resolved. The prior
   "falls back to `tempdir()`" change was insufficient: the default still
   resolved via `system.file("fixtures", package = "TemporalHazard")`, which
   returns the *installed* package directory whenever the package is installed
   (the normal case under `R CMD check`), so the generators still wrote to the
   user library by default. The `output_dir` argument of every fixture
   generator (`.hzr_create_synthetic_golden_fixtures()`,
   `.hzr_create_loglogistic_golden_fixture()`,
   `.hzr_create_lognormal_golden_fixture()`,
   `.hzr_create_multiphase_golden_fixture()`,
   `.hzr_create_c_reference_kul_fixture()`, `.hzr_generate_golden_fixture()`)
   is now **required with no default** — the functions never choose a write
   path themselves. No package code, example, vignette, or test calls these
   generators (they are maintainer-only regeneration helpers); documentation
   shows callers passing `tempdir()` or an explicit project path.

3. **Setting a specific seed within a function.** Now fully resolved.
   `set.seed()` is only ever called inside `if (!is.null(seed))`, where `seed`
   is a user-supplied argument that defaults to `NULL`; the global
   `.Random.seed` is saved beforehand and restored via `on.exit()`. The
   remaining hardcoded `seed = 42` literals in `R/golden_fixtures.R` (stored
   fixture metadata, not `set.seed()` calls) have been replaced with the
   actual `seed` argument, so no specific seed number is set or recorded
   inside any function.

## Test environments

* **Local:** R 4.6.0 on macOS (aarch64-apple-darwin23).
  `R CMD check --as-cran` returns 0 errors, 0 warnings, 0 notes.
* **GitHub Actions matrix:** ubuntu-latest (R-devel / R-release / R-oldrel-1),
  macos-latest (R-release), windows-latest (R-release).
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
