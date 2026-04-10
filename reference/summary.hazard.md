# Summarize a hazard model

Returns a compact summary of a `hazard` object, including model
metadata, fit diagnostics, and coefficient-level statistics when
available.

## Usage

``` r
# S3 method for class 'hazard'
summary(object, ...)
```

## Arguments

- object:

  A `hazard` object.

- ...:

  Unused; for S3 compatibility.

## Value

An object of class `summary.hazard`.

## Examples

``` r
fit <- hazard(time = rexp(30, 0.5), status = rep(1L, 30),
              theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
summary(fit)
#> hazard model summary
#>   observations: 30 
#>   predictors:   0 
#>   dist:         weibull 
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   log-lik:      -45.5527 
#>   evaluations: fn=34, gr=9
#> 
#> Coefficients:
#>     estimate  std_error   z_stat      p_value
#> mu 0.5586338 0.09245866 6.041985 1.522301e-09
#> nu 1.1575846 0.17065867 6.783040 1.176732e-11
```
