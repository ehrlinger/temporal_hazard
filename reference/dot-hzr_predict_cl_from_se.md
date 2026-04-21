# Apply a scale transform + back-transform for CLs

Apply a scale transform + back-transform for CLs

## Usage

``` r
.hzr_predict_cl_from_se(fit, se_nat, level, scale)
```

## Arguments

- fit:

  Point estimate (numeric vector length n).

- se_nat:

  SE on the natural scale (numeric vector length n).

- level:

  Confidence level.

- scale:

  One of "natural", "log", "loglog_survival".

## Value

data.frame with columns fit, se.fit, lower, upper.
