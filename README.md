# hvtiRhazard

`hvtiRhazard` is a native R port of core HAZARD modeling components.

## Development goals

- Preserve numerical behavior where feasible.
- Validate parity against reference HAZARD outputs.
- Keep implementation readable and testable in pure R.

## Development setup

```r
install.packages(c("devtools", "testthat"))
devtools::install_deps(dependencies = TRUE)
devtools::test()
```
