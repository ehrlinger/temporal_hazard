# Select the phase whose log_mu will be solved by conservation

Chooses the phase contributing the largest share of total cumulative
hazard, matching the C HAZARD SETCOE strategy.

## Usage

``` r
.hzr_select_fixmu_phase(
  theta,
  time,
  status,
  phases,
  covariate_counts,
  x_list,
  weights = NULL
)
```

## Arguments

- theta:

  Full parameter vector (internal scale).

- time:

  Numeric vector of follow-up times.

- status:

  Numeric event indicator.

- phases:

  Named list of validated `hzr_phase` objects.

- covariate_counts:

  Named integer vector.

- x_list:

  Named list of per-phase design matrices.

- weights:

  Optional numeric vector of row weights (length n). Defaults to unit
  weights. Applied when summing per-phase cumhaz so that selection
  happens on the same scale as the (weighted) observed event count.

## Value

Character: name of the phase to fix.
