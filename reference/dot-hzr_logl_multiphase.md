# Log-likelihood for multiphase additive hazard model

Log-likelihood for multiphase additive hazard model

## Usage

``` r
.hzr_logl_multiphase(
  theta,
  time,
  status,
  time_lower = NULL,
  time_upper = NULL,
  x = NULL,
  phases,
  covariate_counts,
  x_list,
  return_gradient = FALSE,
  return_hessian = FALSE,
  ...
)
```

## Arguments

- theta:

  Full parameter vector (internal scale).

- time:

  Numeric vector of follow-up times (n).

- status:

  Numeric event indicator: 1 = event, 0 = right-censored, -1 =
  left-censored, 2 = interval-censored.

- time_lower:

  Optional lower bounds for interval censoring.

- time_upper:

  Optional upper bounds for left/interval censoring.

- x:

  Design matrix (unused directly; kept for interface compatibility).

- phases:

  Named list of validated `hzr_phase` objects.

- covariate_counts:

  Named integer vector of per-phase covariate counts.

- x_list:

  Named list of per-phase design matrices.

- return_gradient:

  Logical (ignored; gradient via separate function).

- return_hessian:

  Logical (ignored).

- ...:

  Ignored.

## Value

Scalar log-likelihood. Returns -Inf for infeasible parameters.
