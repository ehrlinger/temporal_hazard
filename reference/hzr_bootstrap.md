# Bootstrap resampling for hazard model coefficients

Resample data with replacement, refit the hazard model on each
replicate, and accumulate coefficient distributions. Returns a tidy data
frame of per-replicate estimates with summary statistics. This is the R
equivalent of the SAS `bootstrap.hazard.sas` macro.

## Usage

``` r
hzr_bootstrap(
  object,
  n_boot = 200L,
  fraction = 1,
  seed = NULL,
  verbose = FALSE
)

# S3 method for class 'hzr_bootstrap'
print(x, digits = 4, ...)
```

## Arguments

- object:

  A fitted `hazard` object (with `fit = TRUE`).

- n_boot:

  Integer: number of bootstrap replicates (default 200).

- fraction:

  Numeric in (0, 1\]: fraction of data to sample per replicate (default
  1.0 for full bootstrap; \< 1 for bagging).

- seed:

  Optional integer random seed for reproducibility.

- verbose:

  Logical; if `TRUE`, print progress every 50 replicates.

- x:

  An `hzr_bootstrap` object.

- digits:

  Number of decimal places for formatting.

- ...:

  Additional arguments (ignored).

## Value

A list with class `"hzr_bootstrap"` containing:

- replicates:

  Data frame with columns `replicate`, `parameter`, and `estimate` – one
  row per parameter per successful replicate.

- summary:

  Data frame with columns `parameter`, `n`, `pct`, `mean`, `sd`, `min`,
  `max`, `ci_lower`, `ci_upper` – one row per parameter.

- n_success:

  Number of successfully converged replicates.

- n_failed:

  Number of replicates that failed to converge.

## See also

[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
for model fitting,
[`vcov.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/vcov.hazard.md)
for Hessian-based standard errors.

## Examples

``` r
# \donttest{
data(avc)
avc <- na.omit(avc)
fit <- hazard(
  survival::Surv(int_dead, dead) ~ age + mal,
  data  = avc,
  dist  = "weibull",
  theta = c(mu = 0.01, nu = 0.5, 0, 0),
  fit   = TRUE
)
bs <- hzr_bootstrap(fit, n_boot = 50, seed = 123)
print(bs)
#> Bootstrap inference for hazard model
#> Replicates: 50 successful, 0 failed
#> 
#>  parameter   n pct   mean     sd     min    max ci_lower ci_upper
#>         mu  50 100 0.0004 0.0004  0.0000 0.0023   0.0000   0.0013
#>         nu  50 100 0.2219 0.0144  0.1963 0.2618   0.1986   0.2462
#>            100 200 0.4201 0.4687 -0.0129 1.3909  -0.0114   1.2998
# }
```
