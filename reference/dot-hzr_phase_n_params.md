# Total number of parameters for a phase

Returns the count of free parameters in the theta sub-vector for one
phase: 1 (log_mu) + n_shape + n_covariates.

## Usage

``` r
.hzr_phase_n_params(phase, n_covariates = 0L)
```

## Arguments

- phase:

  An `hzr_phase` object.

- n_covariates:

  Integer; number of covariate columns this phase uses.

## Value

Integer.
