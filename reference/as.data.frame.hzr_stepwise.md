# Coerce an `hzr_stepwise` result to its selection trace

Returns the `$steps` data frame so downstream tidyverse / data.table
pipelines can work with the trace directly.

## Usage

``` r
# S3 method for class 'hzr_stepwise'
as.data.frame(x, ...)
```

## Arguments

- x:

  An `hzr_stepwise` object.

- ...:

  Ignored.
