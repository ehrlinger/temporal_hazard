# Jacobian of multiphase cumulative-hazard predictions

Only `cumulative_hazard` and `survival` are delegated here – `hazard`
and `linear_predictor` are rejected upstream for multiphase models.

## Usage

``` r
.hzr_predict_jacobian_multiphase(
  theta,
  time,
  phases,
  covariate_counts,
  x_list,
  p
)
```

## Arguments

- theta:

  MLE parameter vector.

- time:

  Prediction times (length n).

- phases:

  Named list of `hzr_phase` objects.

- covariate_counts:

  Named integer vector.

- x_list:

  Named list of per-phase design matrices.

- p:

  Length of theta.

## Value

Numeric n x p Jacobian of H(t\|x).
