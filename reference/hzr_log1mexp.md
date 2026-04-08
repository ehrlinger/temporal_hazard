# Numerically stable log(1 - exp(-x)) for x \> 0

Numerically stable log(1 - exp(-x)) for x \> 0

## Usage

``` r
hzr_log1mexp(x)
```

## Arguments

- x:

  Numeric vector with positive values.

## Value

Numeric vector with element-wise `log(1 - exp(-x))`.

## Examples

``` r
hzr_log1mexp(c(0.01, 0.5, 5))
#> [1] -4.610166019 -0.932752130 -0.006760749
```
