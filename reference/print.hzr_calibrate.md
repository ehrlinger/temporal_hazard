# Print method for hzr_calibrate

Print method for hzr_calibrate

## Usage

``` r
# S3 method for class 'hzr_calibrate'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An `hzr_calibrate` object.

- digits:

  Number of decimal places for formatting.

- ...:

  Additional arguments (ignored).

## Value

The object `x` of class `c("hzr_calibrate", "data.frame")`, invisibly.
The data frame has one row per quantile group and columns: `group`
(group index), `n` (group size), `events` (event count), `mean`, `min`,
`max` (covariate summary within group), `prob` (observed event
probability), `link_value` (transformed probability on the link scale).
When stratified via the `by` argument, a `by` column is also present.
Attributes: `"link"` (the transform applied, e.g. `"logit"`) and
`"groups"` (number of quantile bins).
