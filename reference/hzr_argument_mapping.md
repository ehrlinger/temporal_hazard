# Legacy HAZARD to hvtiRhazard argument mapping

Returns a formal mapping table that defines how legacy SAS
HAZARD/C-style inputs map to `hazard(...)` arguments in this package.

## Usage

``` r
hzr_argument_mapping(include_planned = TRUE)
```

## Arguments

- include_planned:

  Logical; if `FALSE`, only rows marked as implemented are returned.

## Value

A data frame with one row per mapping rule.
