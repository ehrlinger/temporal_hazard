# Extract variance-covariance matrix from hazard model

Returns the estimated variance-covariance matrix of the fitted
coefficients.

## Usage

``` r
# S3 method for class 'hazard'
vcov(object, ...)
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
vcov(fit)
#>              [,1]         [,2]
#> [1,]  0.012518579 -0.005594079
#> [2,] -0.005594079  0.023373248
```
