# Print method for hzr_kaplan

Print method for hzr_kaplan

## Usage

``` r
# S3 method for class 'hzr_kaplan'
print(x, digits = 4, n = 20, ...)
```

## Arguments

- x:

  An `hzr_kaplan` object.

- digits:

  Number of decimal places for formatting.

- n:

  Maximum rows to print (default 20).

- ...:

  Additional arguments (ignored).

## Value

The object `x` of class `c("hzr_kaplan", "data.frame")`, invisibly. The
data frame has one row per event time (or all times when
`event_only = FALSE`) and columns: `time` (follow-up time), `n_risk`
(number at risk), `n_event` (events in interval), `n_censor` (censored
observations in interval), `survival` (Kaplan-Meier survival estimate),
`std_err` (Greenwood standard error on log-hazard scale), `cl_lower`,
`cl_upper` (logit-transformed confidence limits on the survival scale),
`cumhaz` (Nelson-Aalen cumulative hazard), `hazard` (interval hazard
estimate), `density` (estimated event density), `life` (life-table life
expectancy contribution).
