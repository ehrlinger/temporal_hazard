# Extract coefficients from hazard model

Extract coefficients from hazard model

## Usage

``` r
# S3 method for class 'hazard'
coef(object, ...)
```

## Arguments

- object:

  A `hazard` object.

- ...:

  Unused; for S3 compatibility.

## Examples

``` r
fit <- hazard(time = rexp(30, 0.5), status = rep(1L, 30),
              theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
coef(fit)
#> [1]   0.09658355 245.43928875
```
