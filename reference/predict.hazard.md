# Predict from a hazard model object

Produces prediction outputs from a `hazard` object. Supports multiple
prediction types including linear predictor, hazard, survival
probability, and cumulative hazard.

## Usage

``` r
# S3 method for class 'hazard'
predict(
  object,
  newdata = NULL,
  type = c("hazard", "linear_predictor", "survival", "cumulative_hazard"),
  decompose = FALSE,
  ...
)
```

## Arguments

- object:

  A `hazard` object.

- newdata:

  Optional matrix or data frame of predictors. For types requiring time
  (e.g., "survival", "cumulative_hazard"), newdata should include a
  `time` column, or time will be taken from the fitted object's data.

- type:

  Prediction type:

  - `"linear_predictor"`: Linear predictor η = x·β (not available for
    multiphase)

  - `"hazard"`: Hazard scale exp(η) (not available for multiphase)

  - `"survival"`: Survival probability S(t\|x) = exp(-H(t\|x))

  - `"cumulative_hazard"`: Cumulative hazard H(t\|x) at event times

- decompose:

  Logical; if `TRUE` and the model is multiphase, return a data frame
  with per-phase cumulative hazard contributions alongside the total.
  Ignored for single-distribution models. Default `FALSE`.

- ...:

  Unused; included for S3 compatibility.

## Value

Numeric vector of predictions.

## Details

For Weibull models with survival or cumulative_hazard predictions:

- Cumulative hazard: H(t\|x) = (μ·t)^ν · exp(η)

- Survival: S(t\|x) = exp(-H(t\|x))

Time values must be positive and finite. If newdata contains a `time`
column, it will be used; otherwise, the time vector from the fitted
object is used. For models fit with `time_windows`, predictions for
`type = "linear_predictor"` or `"hazard"` also require time values (via
`newdata$time` or fitted-time fallback) so window-specific coefficients
can be selected.

## Examples

``` r
# ── Basic predictions ────────────────────────────────────────────────
set.seed(1)
fit <- hazard(time = rexp(50, 0.3), status = rep(1L, 50),
              theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
predict(fit, type = "survival")
#>  [1] 0.508448836 0.315117497 0.909632142 0.913801355 0.704022509 0.034429522
#>  [7] 0.297904379 0.635876619 0.407732732 0.908685617 0.245861859 0.504733067
#> [13] 0.295096862 0.003729697 0.364948062 0.373066533 0.134498594 0.565324696
#> [19] 0.772692686 0.605268524 0.070992893 0.572924045 0.803140735 0.619329758
#> [25] 0.937239852 0.968075472 0.611315229 0.007471798 0.318195648 0.389681823
#> [31] 0.232969073 0.981596565 0.781843506 0.267479469 0.868327070 0.378412625
#> [37] 0.797695039 0.524953194 0.510431884 0.845615583 0.354514703 0.376046945
#> [43] 0.276615587 0.289749395 0.626388405 0.798022001 0.276332062 0.390676459
#> [49] 0.652266897 0.113525051
predict(fit, newdata = data.frame(time = c(1, 2, 5)),
        type = "cumulative_hazard")
#> [1] 0.2244716 0.5138491 1.5356737

# ── Patient-specific survival curves ─────────────────────────────────
set.seed(1001)
n   <- 180
dat <- data.frame(
  time   = rexp(n, rate = 0.35) + 0.05,
  status = rbinom(n, size = 1, prob = 0.6),
  age    = rnorm(n, mean = 62, sd = 11),
  nyha   = sample(1:4, n, replace = TRUE),
  shock  = rbinom(n, size = 1, prob = 0.18)
)
fit2 <- hazard(
  survival::Surv(time, status) ~ age + nyha + shock,
  data  = dat,
  theta = c(mu = 0.25, nu = 1.10, beta1 = 0, beta2 = 0, beta3 = 0),
  dist  = "weibull", fit = TRUE
)

new_patients <- data.frame(
  time = c(0.5, 1.5, 3.0),
  age  = c(50, 65, 75),
  nyha = c(1, 3, 4),
  shock = c(0, 0, 1)
)
# Compute predictions from the clean covariate frame before adding columns
surv   <- predict(fit2, newdata = new_patients, type = "survival")
cumhaz <- predict(fit2, newdata = new_patients, type = "cumulative_hazard")
new_patients$survival          <- surv
new_patients$cumulative_hazard <- cumhaz
new_patients
#>   time age nyha shock  survival cumulative_hazard
#> 1  0.5  50    1     0 0.9494097        0.05191482
#> 2  1.5  65    3     0 0.7743185        0.25577202
#> 3  3.0  75    4     1 0.5046535        0.68388317

# \donttest{
# ── Grouped survival curves (requires hvtiPlotR) ─────────────────────
if (requireNamespace("hvtiPlotR", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  library(hvtiPlotR)

  t_grid <- seq(0.05, max(dat$time), length.out = 80)
  profiles <- data.frame(
    label = c("Low risk (age 50, NYHA I)",
              "High risk (age 75, NYHA IV)"),
    age   = c(50, 75),
    nyha  = c(1, 4),
    shock = c(0, 1)
  )

  curve_list <- lapply(seq_len(nrow(profiles)), function(i) {
    nd <- data.frame(
      time  = t_grid,
      age   = profiles$age[i],
      nyha  = profiles$nyha[i],
      shock = profiles$shock[i]
    )
    nd$survival <- predict(fit2, newdata = nd, type = "survival") * 100
    nd$profile  <- profiles$label[i]
    nd
  })
  curve_df <- do.call(rbind, curve_list)

  hz_obj <- hv_hazard(
    curve_data   = curve_df,
    x_col        = "time",
    estimate_col = "survival",
    group_col    = "profile"
  )

  plot(hz_obj) +
    ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::labs(
      x     = "Years after surgery",
      y     = "Freedom from death (%)",
      title = "Predicted survival by risk profile"
    ) +
    hv_theme_manuscript()
}

# }
```
