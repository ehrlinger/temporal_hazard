# temporal_hazard Development Roadmap

Progress: M0 complete, M1 in progress.

---

## M0 — Package scaffold (DONE)

- [x] DESCRIPTION, NAMESPACE, LICENSE, NEWS.md, README.md
- [x] roxygen2 documentation workflow
- [x] pkgdown site configuration (`_pkgdown.yml`)
- [x] GitHub Actions: R-CMD-check (multi-platform, multi-version)
- [x] GitHub Actions: test-coverage (covr + Codecov)
- [x] GitHub Actions: pkgdown deploy (GitHub Pages)
- [x] Numerically stable math primitives (`hzr_log1pexp`, `hzr_log1mexp`, `hzr_clamp_prob`)
- [x] Unit tests for math primitives

---

## M1 — API skeleton and migration baseline (IN PROGRESS)

- [x] `hazard()` constructor — accepts time, status, x, theta, dist, control, ...
- [x] `predict.hazard()` S3 method — linear predictor and hazard scale
- [x] `print.hazard()` S3 method
- [x] Formal legacy-to-R argument mapping table (`hzr_argument_mapping()`)
- [x] Unit tests for hazard API and mapping table
- [x] DESCRIPTION notes pure-R-first / Rcpp-later strategy
- [x] SAS-to-R migration vignette (`vignettes/sas-to-r-migration.Rmd`)
  - [x] Side-by-side SAS statement vs R call mapping table
  - [x] Concrete worked example: one end-to-end SAS HAZARD run translated to R
  - [x] Section: Rcpp acceleration roadmap
- [ ] Push repo to GitHub (`git remote add origin` + `git push -u origin main`)
- [ ] Enable GitHub Pages (source: GitHub Actions)

---

## M2 — Model core port (NOT STARTED)

- [ ] Define output object contract for `$fit`:
  - coefficients, log-likelihood, convergence flag, iteration count, gradient norm
- [ ] Port baseline likelihood evaluation for Weibull distribution
- [ ] Port optimizer loop (mirror HAZARD convergence criteria)
- [ ] Expand `predict.hazard()` types: `"survival"`, `"cumulative_hazard"`
- [ ] Parity test harness:
  - [ ] Generate golden fixtures from `hazard` C binary (canonical datasets)
  - [ ] Store in `inst/fixtures/`
  - [ ] `tests/testthat/test-parity-core.R` comparing key metrics with fixed tolerances
- [ ] CI gate: parity tests must pass on Linux before merge
- [ ] Complete R translations in example vignettes (stubs all created in M1):
  - [ ] `vignettes/examples-hz-estimate-hazard.Rmd` (hz.*)
  - [ ] `vignettes/examples-hm-multivariable.Rmd` (hm.*)
  - [ ] `vignettes/examples-hp-prediction.Rmd` (hp.*)
  - [ ] `vignettes/examples-hs-sensitivity.Rmd` (hs.*)
  - [ ] `vignettes/examples-ac-actuarial.Rmd` (ac.*)
  - [ ] `vignettes/examples-bs-bootstrap.Rmd` (bs.*)
  - [ ] `vignettes/examples-lg-exploratory.Rmd` (lg.*)

---

## M3 — Broader model support (NOT STARTED)

- [ ] Additional baseline distributions (log-logistic, log-normal, exponential)
- [ ] Interval censoring support
- [ ] Time-varying coefficients
- [ ] Parity test coverage for edge cases and failure modes
- [ ] Tolerance policy document (per-metric absolute + relative bounds)

---

## M4 — API stabilization and documentation (NOT STARTED)

- [ ] Freeze public API: `hazard()`, `predict.hazard()`, `summary.hazard()`, `coef.hazard()`
- [ ] Add `summary.hazard()` and `coef.hazard()` methods
- [ ] Formula interface: `hazard(Surv(time, status) ~ x1 + x2, data = df)`
- [ ] Full pkgdown reference site with worked examples
- [ ] Migration vignette finalized and reviewed
- [ ] Consider `confint.hazard()` and `vcov.hazard()` if optimizer provides Hessian

---

## M5 — Performance and release candidate (NOT STARTED)

- [ ] Profile against real datasets; identify hot paths
- [ ] Rcpp acceleration for bottlenecks only (likelihood kernel, inner optimizer loop)
  - Interface stays identical; only internal compute changes
- [ ] Add `Rcpp` and `RcppArmadillo` to DESCRIPTION only at this stage
- [ ] CRAN-style submission checks (`--as-cran`)
- [ ] Tag first release `v0.1.0`

---

## Cross-cutting: always required before merge

- [ ] `devtools::test()` passes
- [ ] `R CMD check --as-cran` status: OK (0 warnings, 0 notes)
- [ ] Roxygen outputs current (CI enforces this)
- [ ] Parity tests green on Linux (from M2 onward)
