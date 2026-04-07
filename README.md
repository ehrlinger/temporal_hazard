# temporal_hazard

`temporal_hazard` is a native R implementation of parametric temporal hazard modeling components.

## Current status

- Migration status is tracked in `MIGRATION-STATUS.md`.
- Core modeling is implemented through M3, with major M4 API stabilization now in place.
- Current package validation is clean: tests pass and `R CMD check --no-manual` is OK.

## Development goals

- Preserve numerical behavior where feasible.
- Validate parity against reference HAZARD outputs.
- Keep implementation readable and testable in pure R.

## Documentation entry points

- Start with `vignettes/getting-started.qmd` for a minimal fit-predict workflow.
- Use `vignettes/sas-to-r-migration.qmd` for SAS/HAZARD to R translation guidance.
- Browse the analysis-family examples under `vignettes/` for exploratory, multivariable, prediction, bootstrap, and sensitivity workflows.

## Development setup

```r
install.packages(c("devtools", "roxygen2", "pkgdown", "testthat"))
devtools::install_deps(dependencies = TRUE)
roxygen2::roxygenise()
devtools::test()
pkgdown::build_site()
```

## CI and docs

- GitHub Actions runs multi-platform `R CMD check` on push and pull request.
- CI enforces that roxygen-generated files are current.
- Coverage is published to Codecov.
- pkgdown site is built and deployed from `main` with GitHub Pages.
