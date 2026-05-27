# CRAN submission comments -- TemporalHazard 1.0.3

## Summary

This is a resubmission of 1.0.2 addressing the CRAN reviewer comment that
`R/diagnostics.R` modified `.GlobalEnv`, which is not allowed by CRAN
policy.

1. **Writing to `.GlobalEnv`.** Resolved. In 1.0.2, `hzr_bootstrap()` saved
   and restored the caller's `.Random.seed` via `oldseed <- get0(...)`
   followed by `on.exit({ assign(".Random.seed", oldseed, envir = .GlobalEnv) })`
   (or `rm(...)` when `oldseed` was `NULL`). That hand-rolled save-restore
   wrapper has been removed. When `seed` is supplied the function now
   simply calls `set.seed(seed)` — the documented R API for seeded
   reproducibility — and the `@param seed` documentation notes that the
   caller's RNG state is not restored on exit. The previous reviewer's
   request to avoid setting a *specific* seed inside a function (1.0.1
   feedback) remains satisfied: no `set.seed()` call uses a hardcoded
   value; the only `set.seed()` is conditional on a user-supplied
   argument. No package source file under `R/` references `.GlobalEnv`.

The fixes from 1.0.2 are unchanged: the maintainer-only golden-fixture
generators remain in `data-raw/` (off the CRAN surface);
`.hzr_generate_golden_fixture()` (the one in-package writer) still
requires an explicit `output_dir`; `\value` documentation is present for
every exported function and S3 method.

## Test environments

* **Local:** R 4.6.0 on macOS (aarch64-apple-darwin23).
  `R CMD check --as-cran` returns 0 errors, 0 warnings, 0 notes.
* **GitHub Actions matrix:** ubuntu-latest (R-devel / R-release / R-oldrel-1),
  macos-latest (R-release), windows-latest (R-release).
* **win-builder:** R-devel (`devtools::check_win_devel()`) — 1 NOTE
  (expected; see NOTE disposition below).
* **Reverse-dependency check:** `tools::package_dependencies(reverse = TRUE)`
  returns 0.

## NOTE disposition

Local `R CMD check --as-cran` is clean (0 errors, 0 warnings, 0 notes).
CRAN's incoming-feasibility check (and win-builder) report a single
expected NOTE — the local check does not exercise the incoming-feasibility
/ aspell step, which is why it is surfaced only by CRAN/win-builder:

* **New submission.** TemporalHazard is not yet on CRAN. *(Remove this
  bullet once the package is accepted — subsequent releases are updates,
  not new submissions.)*
* **Possibly misspelled words in DESCRIPTION:** `Naftel`, `Rajeswaran`
  (cited authors, appearing with their `<doi:...>` references), `UAB`
  (acronym — University of Alabama at Birmingham), `et`, `al`
  ("et al." citation), and `multiphase` (the core statistical term this
  package implements). None are misspellings; all are intentional proper
  nouns, an acronym, a standard citation, and a domain term. No action
  taken.

(`inst/WORDLIST` already lists these words for the `spelling` package's
unit test, but the CRAN incoming aspell check has no package-level
suppression mechanism, so this NOTE is informational and recurs by design
on every submission. Keep this section in sync with the actual
CRAN/win-builder output, not the local `--as-cran` result.)

## Background

* Package contains five clinical reference datasets (~120 KB compressed) used
  in vignettes and parity tests against the original C/SAS HAZARD program.
* Vignettes use Quarto (`VignetteBuilder: quarto`).
* A subset of tests (SAS parity, multi-start multiphase fits, OMC raw-file
  derivation) are guarded with `skip_on_cran()` to stay within CRAN runtime
  budgets; the full suite runs on GitHub Actions.
* DOI URLs in vignettes (doi.org) may return HTTP 403 to automated checkers
  but resolve correctly in browsers.
