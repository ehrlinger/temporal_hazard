# CRAN submission comments — TemporalHazard 0.9.1

## Summary

First submission. TemporalHazard provides native R implementations of the
multiphase parametric hazard model (Blackstone, Naftel, Turner 1986) for
temporal decomposition of survival data into additive hazard phases.

## Test environments

* local macOS (aarch64-apple-darwin), R 4.x.x
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1),
  macos-latest (release), windows-latest (release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes

* This is a new submission.
* Package contains five clinical reference datasets (total ~120 KB compressed)
  used in vignettes and parity tests against the original C/SAS HAZARD program.
* Vignettes use Quarto (`VignetteBuilder: quarto`).
* DOI URLs in vignettes (doi.org) may return HTTP 403 to automated checkers
  but resolve correctly in browsers.
