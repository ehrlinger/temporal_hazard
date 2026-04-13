# Transform multiphase theta from internal to user-facing scale

Converts log_mu -\> mu, log_t_half -\> t_half; nu, m, beta pass through.

## Usage

``` r
.hzr_multiphase_theta_user(theta, phases, covariate_counts)
```

## Arguments

- theta:

  Full parameter vector (internal scale).

- phases:

  Named list of validated `hzr_phase` objects.

- covariate_counts:

  Named integer vector of per-phase covariate counts.

## Value

Named numeric vector on the reporting scale.
