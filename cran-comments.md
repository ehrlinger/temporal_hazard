# CRAN submission comments -- TemporalHazard 1.0.0

## Summary

First submission. TemporalHazard provides native R implementations of the
multiphase parametric hazard model (Blackstone, Naftel, Turner 1986) for
temporal decomposition of survival data into additive hazard phases.

Version 1.0.0 is the production release. Key capabilities:

* Multiphase additive hazard fitting (`hazard()`) with CDF, constant, and
  generalized three-parameter (G3) phase types
* Observation weights, start-stop (`Surv(start, stop, event)`) left-truncated
  and repeating-event data, Conservation of Events theorem
* Seven SAS-macro equivalent diagnostic/utility functions: `hzr_kaplan`,
  `hzr_nelson`, `hzr_gof`, `hzr_deciles`, `hzr_calibrate`, `hzr_bootstrap`,
  `hzr_competing_risks`
* Delta-method confidence limits on all prediction types
* SAS parity test suite (54 expectations / 6 model fits) validating numerical
  agreement with the original C/SAS HAZARD program across four clinical datasets

## Test environments

* local macOS (aarch64-apple-darwin), R 4.6.x
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1),
  macos-latest (release), windows-latest (release)

## R CMD check results

0 errors | 0 warnings | 1 note

The note is:

```
checking CRAN incoming feasibility ...
  Maintainer: 'John Ehrlinger <john.ehrlinger@gmail.com>'

New submission
```

This is expected for a new submission.

## Notes

* This is a new submission.
* Package contains five clinical reference datasets (total ~120 KB compressed)
  used in vignettes and parity tests against the original C/SAS HAZARD program.
* Vignettes use Quarto (`VignetteBuilder: quarto`).
* A subset of tests (SAS parity, multi-start multiphase fits, OMC raw-file
  derivation) are guarded with `skip_on_cran()` to stay within CRAN runtime
  budgets; the full suite runs on GitHub Actions.
* DOI URLs in vignettes (doi.org) may return HTTP 403 to automated checkers
  but resolve correctly in browsers.
  