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

## See also

[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
for model fitting,
[`summary.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/summary.hazard.md)
for model summaries,
[`hzr_phase()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase.md)
for multiphase temporal shapes.

[`vignette("prediction-visualization")`](https://ehrlinger.github.io/temporal_hazard/articles/prediction-visualization.md)
for detailed prediction workflows including decomposed hazard plots and
patient-specific curves.

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
# ── Grouped survival curves ───────────────────────────────────────
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

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

  ggplot(curve_df, aes(time, survival, colour = profile)) +
    geom_line() +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "Months after surgery",
         y = "Freedom from death (%)",
         title = "Predicted survival by risk profile",
         colour = NULL) +
    theme_minimal()
}

# }

# \donttest{
# ── Multiphase predictions with decomposition ────────────────────
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

t_grid <- seq(0.01, max(dat$time) * 0.9, length.out = 100)
nd     <- data.frame(time = t_grid)

# Overall survival
predict(fit_mp, newdata = nd, type = "survival")
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.99963159   0.97022789   0.93193549   0.89108549   0.84973756   0.80894020 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.76928444   0.73110930   0.69460102   0.65984848   0.62687690   0.59566952 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.56618210   0.53835279   0.51210901   0.48737225   0.46406139   0.44209506 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.42139326   0.40187840   0.38347603   0.36611526   0.34972900   0.33425403 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.31963102   0.30580442   0.29272234   0.28033640   0.26860153   0.25747584 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.24692035   0.23689886   0.22737774   0.21832577   0.20971395   0.20151535 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.19370496   0.18625955   0.17915752   0.17237880   0.16590473   0.15971792 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.15380223   0.14814257   0.14272492   0.13753619   0.13256417   0.12779745 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.12322538   0.11883799   0.11462597   0.11058058   0.10669366   0.10295753 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.09936500   0.09590930   0.09258410   0.08938341   0.08630162   0.08333343 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.08047384   0.07771813   0.07506188   0.07250085   0.07003109   0.06764883 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.06535049   0.06313271   0.06099227   0.05892612   0.05693138   0.05500530 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.05314524   0.05134873   0.04961338   0.04793693   0.04631722   0.04475217 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.04323981   0.04177827   0.04036572   0.03900044   0.03768078   0.03640514 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.03517200   0.03397989   0.03282742   0.03171323   0.03063603   0.02959457 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.02858765   0.02761411   0.02667285   0.02576278   0.02488289   0.02403216 
#> early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.02320965   0.02241442   0.02164559   0.02090229 

# Per-phase decomposed cumulative hazard
decomp <- predict(fit_mp, newdata = nd,
                  type = "cumulative_hazard", decompose = TRUE)
head(decomp)
#>        time        total        early         late
#> 1 0.0100000 0.0003684827 0.0003042602 6.422247e-05
#> 2 0.3177112 0.0302242973 0.0238913236 6.332974e-03
#> 3 0.6254224 0.0704916794 0.0549308704 1.556081e-02
#> 4 0.9331336 0.1153149089 0.0888509069 2.646400e-02
#> 5 1.2408448 0.1628277280 0.1241998078 3.862792e-02
#> 6 1.5485560 0.2120302841 0.1602039971 5.182629e-02
# }
```
