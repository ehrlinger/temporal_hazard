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
#>   log-lik:      31843.7 
#>   evaluations: fn=126, gr=126
#>   message:      ERROR: ABNORMAL_TERMINATION_IN_LNSRCH 
#> 
#> Coefficients:
#>           estimate std_error z_stat p_value
#> mu      0.02629409        NA     NA      NA
#> nu    194.92997141        NA     NA      NA
#> beta1   2.25357255        NA     NA      NA
#> beta2   0.49084734        NA     NA      NA
#> beta3   0.02367525        NA     NA      NA
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
#> 1  0.5  50    1     0         113.1695      1.408895e+49        1
#> 2  1.5  65    3     0         147.9548      1.802748e+64        1
#> 3  3.0  75    4     1         171.0050      1.847271e+74        1
#>   cumulative_hazard
#> 1      0.000000e+00
#> 2     3.667761e-210
#> 3     1.797898e-141
```

## Numerical helpers

The numerical helper functions remain available directly when you need
stable log-scale calculations for custom work or debugging.

``` r
TemporalHazard::hzr_log1pexp(c(-2, 0, 2))
#> [1] 0.1269280 0.6931472 2.1269280
```
