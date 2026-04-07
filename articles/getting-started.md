# Getting Started with TemporalHazard

This vignette shows the minimal workflow for fitting a parametric hazard
model, summarizing the fit, and generating predictions.

``` r
if (!requireNamespace("TemporalHazard", quietly = TRUE)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("TemporalHazard is unavailable and pkgload is not installed.", call. = FALSE)
  }
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
}
library(survival)

set.seed(1001)
n <- 180
dat <- data.frame(
  time = rexp(n, rate = 0.35) + 0.05,
  status = rbinom(n, size = 1, prob = 0.6),
  age = rnorm(n, mean = 62, sd = 11),
  nyha = sample(1:4, n, replace = TRUE),
  shock = rbinom(n, size = 1, prob = 0.18)
)

fit <- TemporalHazard::hazard(
  Surv(time, status) ~ age + nyha + shock,
  data = dat,
  theta = c(mu = 0.25, nu = 1.10, beta1 = 0, beta2 = 0, beta3 = 0),
  dist = "weibull",
  fit = TRUE,
  control = list(maxit = 300)
)

summary(fit)
#> hazard model summary
#>   observations: 180 
#>   predictors:   3 
#>   dist:         weibull 
#>   engine:       native-r-m2 
#>   converged:    FALSE 
#>   log-lik:      28299.2 
#>   evaluations: fn=112, gr=112
#>   message:      ERROR: ABNORMAL_TERMINATION_IN_LNSRCH 
#> 
#> Coefficients:
#>           estimate std_error z_stat p_value
#> mu      0.02056897        NA     NA      NA
#> nu    209.33313265        NA     NA      NA
#> beta1   1.55485771        NA     NA      NA
#> beta2   0.49007717        NA     NA      NA
#> beta3   0.02430348        NA     NA      NA
```

## Prediction workflow

``` r
new_patients <- data.frame(
  time = c(0.5, 1.5, 3.0),
  age = c(50, 65, 75),
  nyha = c(1, 3, 4),
  shock = c(0, 0, 1)
)

pred_input <- new_patients

new_patients$linear_predictor <- predict(fit, newdata = pred_input, type = "linear_predictor")
new_patients$hazard_multiplier <- predict(fit, newdata = pred_input, type = "hazard")
new_patients$survival <- predict(fit, newdata = pred_input, type = "survival")
new_patients$cumulative_hazard <- predict(fit, newdata = pred_input, type = "cumulative_hazard")

new_patients
#>   time age nyha shock linear_predictor hazard_multiplier survival
#> 1  0.5  50    1     0         78.23296      9.465509e+33        1
#> 2  1.5  65    3     0        102.53598      3.394779e+44        1
#> 3  3.0  75    4     1        118.59894      3.212664e+51        1
#>   cumulative_hazard
#> 1      0.000000e+00
#> 2     1.959013e-272
#> 3     1.921513e-202
```

## Numerical helpers

The numerical helper functions remain available directly when you need
stable log-scale calculations for custom work or debugging.

``` r
TemporalHazard::hzr_log1pexp(c(-2, 0, 2))
#> [1] 0.1269280 0.6931472 2.1269280
```
