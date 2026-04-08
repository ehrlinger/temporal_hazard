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
#>   log-lik:      4350.07 
#>   evaluations: fn=167, gr=167
#>   message:      CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH 
#> 
#> Coefficients:
#>      estimate    std_error  z_stat p_value
#> mu   0.102828 1.287868e-08 7984361       0
#> nu 312.346556 3.911979e-05 7984361       0
```
