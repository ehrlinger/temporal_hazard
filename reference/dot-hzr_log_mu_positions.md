# Find the position of each phase's log_mu in the full theta vector

Find the position of each phase's log_mu in the full theta vector

## Usage

``` r
.hzr_log_mu_positions(phases, covariate_counts)
```

## Arguments

- phases:

  Named list of validated `hzr_phase` objects.

- covariate_counts:

  Named integer vector of per-phase covariate counts.

## Value

Named integer vector: position of log_mu for each phase.
