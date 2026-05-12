# Extract the captured console trace from an `hzr_stepwise` fit

Every run of
[`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
records the header, per-step lines, and final summary regardless of the
`trace` flag. This accessor returns the full character vector for
display or logging.

## Usage

``` r
stepwise_trace(fit)
```

## Arguments

- fit:

  An `hzr_stepwise` object.

## Value

Character vector, one element per console line.

## Examples

``` r
data(avc)
avc <- na.omit(avc)
base <- hazard(survival::Surv(int_dead, dead) ~ age,
               data = avc, dist = "weibull", fit = TRUE)
# \donttest{
sw <- hzr_stepwise(base, scope = ~ age + nyha,
                   data = avc, direction = "forward",
                   control = list(n_starts = 1))
#> Stepwise selection (direction = forward, criterion = wald, slentry = 0.30, slstay = 0.20)
#> 
#> Warning: Stepwise forward: candidate refit failed for nyha.
#> (no further action after 0 steps)
#> 
#> Final model: 0 covariates, logLik = NA, AIC = NA
cat(stepwise_trace(sw), sep = "\n")
#> Stepwise selection (direction = forward, criterion = wald, slentry = 0.30, slstay = 0.20)
#> 
#> (no further action after 0 steps)
#> 
#> Final model: 0 covariates, logLik = NA, AIC = NA
# }
```
