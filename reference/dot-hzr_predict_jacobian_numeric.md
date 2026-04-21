# Numeric Jacobian of a single-distribution prediction w.r.t. theta

Numeric Jacobian of a single-distribution prediction w.r.t. theta

## Usage

``` r
.hzr_predict_jacobian_numeric(predict_fn, theta)
```

## Arguments

- predict_fn:

  Function(theta) -\> numeric vector of length n.

- theta:

  MLE parameter vector.

## Value

Numeric n x p Jacobian.
