# Parse Surv() formula for hazard modeling

Extracts time, status, time_lower, time_upper, and predictors from a
formula of the form Surv(time, status) ~ x1 + x2 + ... Supports
right-censored, left-censored, and interval-censored data.

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
