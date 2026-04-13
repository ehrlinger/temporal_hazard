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
#>   0.99945805   0.96587821   0.92523393   0.88291701   0.84056173   0.79898572 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.75864953   0.71981990   0.68264609   0.64720169   0.61350985   0.58155943 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.55131597   0.52272923   0.49573860   0.47027700   0.44627364   0.42365616 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.40235202   0.38228970   0.36339937   0.34561348   0.32886711   0.31309819 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.29824764   0.28425940   0.27108046   0.25866078   0.24695324   0.23591350 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.22549994   0.21567352   0.20639764   0.19763803   0.18936259   0.18154132 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.17414614   0.16715081   0.16053080   0.15426317   0.14832651   0.14270076 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.13736723   0.13230842   0.12750797   0.12295059   0.11862198   0.11450877 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.11059845   0.10687932   0.10334041   0.09997148   0.09676290   0.09370569 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.09079140   0.08801213   0.08536046   0.08282944   0.08041253   0.07810361 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.07589693   0.07378706   0.07176894   0.06983778   0.06798908   0.06621861 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.06452240   0.06289668   0.06133792   0.05984280   0.05840816   0.05703104 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.05570863   0.05443829   0.05321750   0.05204390   0.05091524   0.04982936 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.04878425   0.04777795   0.04680860   0.04587441   0.04497366   0.04410466 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.04326575   0.04245531   0.04167168   0.04091320   0.04017815   0.03946472 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.03877098   0.03809482   0.03743392   0.03678567   0.03614711   0.03551480 
#> early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.03488478   0.03425242   0.03361226   0.03295793 

# Per-phase decomposed cumulative hazard
decomp <- predict(fit_mp, newdata = nd,
                  type = "cumulative_hazard", decompose = TRUE)
head(decomp)
#>        time        total        early late
#> 1 0.0100000 0.0005420998 0.0005420998    0
#> 2 0.3177112 0.0347175326 0.0347175326    0
#> 3 0.6254224 0.0777086758 0.0777086758    0
#> 4 0.9331336 0.1245240653 0.1245240653    0
#> 5 1.2408448 0.1736848788 0.1736848788    0
#> 6 1.5485560 0.2244122111 0.2244122111    0
# }
```
