# Compute total cumulative hazard H(t\|x) for multiphase model

Compute total cumulative hazard H(t\|x) for multiphase model

## Usage

``` r
.hzr_multiphase_cumhaz(
  time,
  theta,
  phases,
  covariate_counts,
  x_list,
  per_phase = FALSE
)
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

- per_phase:

  Logical; if TRUE return a list with total and per-phase contributions.

## Value

If `per_phase = FALSE`: numeric vector of length n (total H(t\|x)). If
`per_phase = TRUE`: named list with `$total` and one element per phase.
