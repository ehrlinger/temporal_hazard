# Per-row standard errors via the delta method

Per-row standard errors via the delta method

## Usage

``` r
.hzr_predict_se_from_jacobian(J, vcov)
```

## Arguments

- J:

  Jacobian (n x p).

- vcov:

  p x p variance-covariance matrix.

## Value

Numeric vector of length n.
