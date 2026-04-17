# Print method for hazard summary objects

Formatted console display of
[`summary.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/summary.hazard.md)
output: distribution, phase list (for multiphase), coefficient table
with standard errors, and log-likelihood. S3 dispatch only – users call
`print(summary(fit))` rather than invoking this directly.

## Usage

``` r
# S3 method for class 'summary.hazard'
print(x, ...)
```

## Arguments

- x:

  A `summary.hazard` object returned by
  [`summary.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/summary.hazard.md).

- ...:

  Additional arguments (ignored).

## Value

`x`, invisibly.
