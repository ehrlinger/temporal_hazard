# Test if an object is an hzr_phase

Test if an object is an hzr_phase

## Usage

``` r
is_hzr_phase(x)
```

## Arguments

- x:

  Object to test.

## Value

Logical scalar.

## Examples

``` r
is_hzr_phase(hzr_phase("cdf"))
#> [1] TRUE
is_hzr_phase("not a phase")
#> [1] FALSE
```
