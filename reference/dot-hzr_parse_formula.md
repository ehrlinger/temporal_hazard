# Parse Surv() formula for hazard modeling

Extracts time, status, time_lower, time_upper, and predictors from a
formula of the form `Surv(time, status) ~ x1 + x2 + ...`. Supports
right-censored, left-censored, interval-censored, and counting-process
(start-stop) data.

## Usage

``` r
.hzr_parse_formula(formula, data)
```

## Arguments

- formula:

  A formula object with Surv() on the LHS.

- data:

  A data frame containing variables referenced in the formula.

## Value

A list with elements: time, status, time_lower, time_upper, x

## Details

For counting-process (start-stop) data, use `Surv(start, stop, event)`.
The start times are returned as `time_lower` and stop times as `time`,
enabling the likelihood to compute `H(stop) - H(start)` per epoch.
