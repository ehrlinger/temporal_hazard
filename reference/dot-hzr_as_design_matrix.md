# Coerce x to a validated numeric design matrix

Accepts data.frame or numeric matrix input; validates dimensions and
finiteness. Called by hazard() to normalise the x argument before any
downstream use.

## Usage

``` r
.hzr_as_design_matrix(x, n = NULL)
```

## Arguments

- x:

  A numeric matrix or data frame with numeric columns.

- n:

  Expected row count (length of time/status vectors); checked if
  non-NULL.

## Value

A numeric matrix with column names preserved (or added if absent).
