# hvtiRhazard

`hvtiRhazard` is a native R port of core HAZARD modeling components.

## Development goals

- Preserve numerical behavior where feasible.
- Validate parity against reference HAZARD outputs.
- Keep implementation readable and testable in pure R.

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
