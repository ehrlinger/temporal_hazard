# Decile-of-risk calibration

Partition observations into groups (default 10) by predicted risk and
compare observed vs. expected event counts in each group. Good
calibration means the two track each other across the risk spectrum.

## Usage

``` r
hzr_deciles(object, time, groups = 10L, status = NULL, event_time = NULL)
```

## Arguments

- object:

  A fitted `hazard` object (with `fit = TRUE`).

- time:

  Numeric scalar: the time point at which to evaluate predicted survival
  / cumulative hazard. For example, `time = 12` evaluates 12-month
  predictions.

- groups:

  Integer: number of risk groups (default 10 for deciles).

- status:

  Optional numeric vector of event indicators (1 = event, 0 = censored).
  If `NULL` (default), extracted from the fitted object's stored data.

- event_time:

  Optional numeric vector of observed event/censoring times. If `NULL`
  (default), extracted from the fitted object.

## Value

A data frame with one row per risk group and columns:

- group:

  Integer group label (1 = lowest risk).

- n:

  Number of observations in the group.

- events:

  Observed event count by the requested horizon (event indicator = 1 and
  event time \<= `time`).

- expected:

  Expected event count by the requested horizon, computed as the sum of
  individual event probabilities (\\1 - S(time)\\).

- observed_rate:

  Observed event rate (events / n).

- expected_rate:

  Expected event rate (expected / n).

- chi_sq:

  Chi-square contribution: (events - expected)^2 / expected.

- p_value:

  Upper-tail p-value from the chi-square test for this group (1 df).

- mean_survival:

  Mean predicted survival probability in the group.

- mean_cumhaz:

  Mean predicted cumulative hazard in the group.

An attribute `"overall"` is attached with the overall chi-square
statistic, degrees of freedom, and p-value.

## Details

This implements the workflow of the SAS `deciles.hazard.sas` macro.
Patients are ranked by predicted cumulative hazard at a specified time
point, grouped into quantile bins, and each bin is tested with a
chi-square goodness-of-fit statistic. Subjects censored before the
requested horizon are excluded from the observed-vs-expected comparison.

## See also

[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
for the prediction types used internally.

## Examples

``` r
# \donttest{
data(avc)
avc <- na.omit(avc)
fit <- hazard(
  survival::Surv(int_dead, dead) ~ age + mal,
  data  = avc,
  dist  = "weibull",
  theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
  fit   = TRUE
)
cal <- hzr_deciles(fit, time = 120)
print(cal)
#> Decile-of-risk calibration at time = 120 
#> Included 86 observations (excluded 219 censored before horizon).
#> 10 groups, 66 observed events, 29.4 expected
#> 
#>  group n events expected observed_rate expected_rate chi_sq  p_value
#>      1 8      3    0.441         0.375        0.0552 14.800 0.000118
#>      2 9      3    1.450         0.333        0.1610  1.670 0.196000
#>      3 8      3    1.780         0.375        0.2220  0.841 0.359000
#>      4 9      8    2.640         0.889        0.2940 10.900 0.000984
#>      5 9      9    2.930         1.000        0.3260 12.500 0.000399
#>      6 8      8    2.690         1.000        0.3360 10.500 0.001190
#>      7 9      9    3.120         1.000        0.3460 11.100 0.000857
#>      8 8      7    3.690         0.875        0.4610  2.980 0.084200
#>      9 9      7    5.210         0.778        0.5780  0.618 0.432000
#>     10 9      9    5.460         1.000        0.6070  2.300 0.130000
#>  mean_survival mean_cumhaz
#>          0.945      0.0579
#>          0.839      0.1760
#>          0.778      0.2510
#>          0.706      0.3480
#>          0.674      0.3950
#>          0.664      0.4090
#>          0.654      0.4250
#>          0.539      0.6240
#>          0.422      0.8640
#>          0.393      0.9330
#> 
#> Overall: chi-sq = 68.3 on 9 df, p = 3.34e-11 
# }
```
