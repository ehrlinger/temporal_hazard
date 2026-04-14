# Summarize a hazard model

Returns a compact summary of a `hazard` object, including model
metadata, fit diagnostics, and coefficient-level statistics when
available.

## Usage

``` r
# S3 method for class 'hazard'
summary(object, ...)
```

## Arguments

- object:

  A `hazard` object.

- ...:

  Unused; for S3 compatibility.

## Value

An object of class `summary.hazard`.

## See also

[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
for model fitting,
[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
for predictions.

[`vignette("fitting-hazard-models")`](https://ehrlinger.github.io/temporal_hazard/articles/fitting-hazard-models.md)
for fitting workflows,
[`vignette("inference-diagnostics")`](https://ehrlinger.github.io/temporal_hazard/articles/inference-diagnostics.md)
for bootstrap CIs and diagnostics.

## Examples

``` r
# ── Single-phase Weibull summary ────────────────────────────────────
fit <- hazard(time = rexp(30, 0.5), status = rep(1L, 30),
              theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
summary(fit)
#> hazard model summary
#>   observations: 30 
#>   predictors:   0 
#>   dist:         weibull 
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   log-lik:      -57.588 
#>   evaluations: fn=17, gr=5
#> 
#> Coefficients:
#>     estimate  std_error   z_stat      p_value
#> mu 0.4017570 0.07799612 5.150987 2.591190e-07
#> nu 0.9772992 0.14742788 6.628998 3.379725e-11

# \donttest{
# ── Multiphase model summary ────────────────────────────────────────
set.seed(42)
n   <- 200
dat <- data.frame(
  time   = rexp(n, rate = 0.25) + 0.01,
  status = rbinom(n, size = 1, prob = 0.65)
)
fit_mp <- hazard(
  survival::Surv(time, status) ~ 1,
  data   = dat,
  dist   = "multiphase",
  phases = list(
    early = hzr_phase("cdf",      t_half = 0.5, nu = 2, m = 0),
    late  = hzr_phase("cdf",      t_half = 5,   nu = 1, m = 0)
  ),
  fit     = TRUE,
  control = list(n_starts = 3, maxit = 500)
)
summary(fit_mp)
#> Warning: NaNs produced
#> Multiphase hazard model (2 phases)
#>   observations: 200 
#>   predictors:   0 
#>   dist:         multiphase 
#>   phase 1:      early - cdf (early risk)
#>   phase 2:      late - cdf (late risk)
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   log-lik:      -393.584 
#>   evaluations: fn=88, gr=21
#> 
#> Coefficients (internal scale):
#> 
#>   Phase: early (cdf)
#>                  estimate std_error    z_stat     p_value
#>   log_mu      0.239134485       NaN        NA          NA
#>   log_t_half  1.818478590       NaN        NA          NA
#>   nu          0.004353297       NaN        NA          NA
#>   m          -0.787666114 0.1325287 -5.943362 2.79235e-09
#> 
#>   Phase: late (cdf)
#>                 estimate std_error     z_stat    p_value
#>   log_mu      3.74165145 5.9132228  0.6327601 0.52689032
#>   log_t_half  2.39068098 5.9151369  0.4041633 0.68609267
#>   nu         -0.75329358 0.3993425 -1.8863346 0.05924987
#>   m           0.02353084 0.1145096  0.2054924 0.83718746
# }
```
