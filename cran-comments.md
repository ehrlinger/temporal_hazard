# CRAN submission comments -- TemporalHazard 0.9.4

## Summary

First submission. TemporalHazard provides native R implementations of the
multiphase parametric hazard model (Blackstone, Naftel, Turner 1986) for
temporal decomposition of survival data into additive hazard phases.

Since the 0.9.1 draft, 0.9.4 adds observation weights (`WEIGHT` statement
equivalent), repeating events via `Surv(start, stop, event)` start-stop
notation, and the Conservation of Events theorem integrated into the
multiphase optimizer. Seven SAS-macro equivalent diagnostic/utility
functions (`hzr_kaplan`, `hzr_nelson`, `hzr_gof`, `hzr_deciles`,
`hzr_calibrate`, `hzr_bootstrap`, `hzr_competing_risks`) are exported.

## Test environments

* local macOS (aarch64-apple-darwin), R 4.5.x
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1),
  macos-latest (release), windows-latest (release)

## R CMD check results

0 errors | 0 warnings | 1 note (new submission)

## Notes

* This is a new submission.
* Package contains five clinical reference datasets (total ~120 KB compressed)
  used in vignettes and parity tests against the original C/SAS HAZARD program.
* Vignettes use Quarto (`VignetteBuilder: quarto`).
* A subset of tests (C-binary parity, multi-start multiphase parity against
  pre-computed reference fixtures) are guarded with `skip_on_cran()` to stay
  within CRAN runtime budgets; the full suite runs on GitHub Actions.
* DOI URLs in vignettes (doi.org) may return HTTP 403 to automated checkers
  but resolve correctly in browsers.
