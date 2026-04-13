# Validate a list of phase specifications

Checks that `phases` is a non-empty named list of `hzr_phase` objects.
Auto-names unnamed phases as `phase_1`, `phase_2`, etc.

## Usage

``` r
.hzr_validate_phases(phases)
```

## Arguments

- phases:

  A list of `hzr_phase` objects.

## Value

The validated (and possibly auto-named) list, invisibly.
