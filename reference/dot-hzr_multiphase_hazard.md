# Compute total instantaneous hazard h(t\|x) for multiphase model

Compute total instantaneous hazard h(t\|x) for multiphase model

## Usage

``` r
.hzr_multiphase_hazard(time, theta, phases, covariate_counts, x_list)
```

## Arguments

- time:

  Numeric vector of times (n).

- theta:

  Full parameter vector (internal scale).

- phases:

  Named list of validated `hzr_phase` objects.

- covariate_counts:

  Named integer vector of covariate counts per phase.

- x_list:

  Named list of design matrices, one per phase. Each element is an n x
  p_j matrix or NULL.

## Value

Numeric vector of length n.
