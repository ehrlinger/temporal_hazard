# Build a logical mask of free (non-fixed) parameters in the full theta vector

Returns a logical vector the same length as the full theta where `TRUE`
means the parameter is free (to be optimized) and `FALSE` means it is
held fixed at its starting value.

## Usage

``` r
.hzr_phase_free_mask(phases, covariate_counts)
```

## Arguments

- phases:

  Named list of validated `hzr_phase` objects.

- covariate_counts:

  Named integer vector of per-phase covariate counts.

## Value

Logical vector of length `sum(.hzr_phase_n_params(...))`.
