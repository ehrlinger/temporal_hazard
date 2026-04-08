# Numerically stable log(1 + exp(x))

Numerically stable log(1 + exp(x))

## Usage

``` r
hzr_log1pexp(x)
```

## Arguments

- x:

  Numeric vector.

## Value

Numeric vector with element-wise `log(1 + exp(x))`.

## Examples

``` r
hzr_log1pexp(c(-50, 0, 50))
#> [1] 1.928750e-22 6.931472e-01 5.000000e+01
```
