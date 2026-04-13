# TemporalHazard

<!-- badges: start -->
[![version](https://img.shields.io/badge/version-0.9.0-blue.svg)](https://github.com/ehrlinger/temporal_hazard)
[![R-CMD-check](https://github.com/ehrlinger/temporal_hazard/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ehrlinger/temporal_hazard/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ehrlinger/temporal_hazard/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ehrlinger/temporal_hazard?branch=main)
[![pkgdown site](https://img.shields.io/badge/docs-pkgdown-blue)](https://ehrlinger.github.io/temporal_hazard/)
![active](https://www.repostatus.org/badges/latest/active.svg)
<!-- badges: end -->

**TemporalHazard** is a pure-R implementation of the multiphase parametric
hazard model of Blackstone, Naftel, and Turner (1986). It decomposes the
overall hazard of an event into additive temporal phases --- early, constant,
and late --- each governed by the generalized temporal decomposition family.
This structure captures real clinical risk patterns that standard single-distribution
models (Weibull, log-normal) cannot represent.

## Why multiphase?

After cardiac surgery, the risk of death is not constant. It starts high in
the immediate post-operative period (early phase), settles to a low background
rate (constant phase), and eventually rises again as patients age (late phase).
A single Weibull curve forces a monotone shape; a multiphase model captures all
three regimes simultaneously.

<img src="man/figures/readme-hazard-phases.png" width="100%" alt="Additive phase decomposition showing early, constant, and late hazard components summing to the total hazard curve" />

The resulting survival curve closely tracks the nonparametric Kaplan-Meier estimate
while providing a smooth, parametric representation that supports covariate
adjustment, prediction, and extrapolation.

<img src="man/figures/readme-survival.png" width="100%" alt="Multiphase parametric survival curve overlaid on the Kaplan-Meier estimate from the AVC dataset" />

## Quick start

```r
install.packages("remotes")
remotes::install_github("ehrlinger/temporal_hazard")
```

### Single-phase model

```r
library(TemporalHazard)

fit <- hazard(
  survival::Surv(time, status) ~ age + nyha,
  data  = dat,
  dist  = "weibull",
  theta = c(mu = 0.25, nu = 1.1, beta1 = 0, beta2 = 0),
  fit   = TRUE
)
summary(fit)
predict(fit, newdata = new_patients, type = "survival")
```

### Multiphase model

```r
fit_mp <- hazard(
  survival::Surv(int_dead, dead) ~ 1,
  data   = avc,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf",      t_half = 0.5, nu = 1, m = 1),
    constant = hzr_phase("constant"),
    late     = hzr_phase("cdf",      t_half = 10,  nu = 1, m = 1)
  ),
  fit = TRUE
)
summary(fit_mp)

# Per-phase decomposition of cumulative hazard
predict(fit_mp, newdata = grid, type = "cumulative_hazard", decompose = TRUE)
```

Each phase is specified with `hzr_phase()`, which sets the temporal shape
type and starting values. The optimizer estimates both the phase-specific
scale parameters and shape parameters jointly.

## Documentation

- **[Getting Started](https://ehrlinger.github.io/temporal_hazard/articles/getting-started.html)** --- single-phase and multiphase fit-predict workflows with visualizations.
- **[Mathematical Foundations](https://ehrlinger.github.io/temporal_hazard/articles/mf-mathematical-foundations.html)** --- the generalized decomposition, additive hazard model, censoring likelihood, and time-varying covariates.
- **[Package Architecture](https://ehrlinger.github.io/temporal_hazard/articles/ar-architecture.html)** --- internal design, golden fixtures, and dataset catalog.
- **[SAS-to-R Migration](https://ehrlinger.github.io/temporal_hazard/articles/sas-to-r-migration.html)** --- statement-by-statement mapping from SAS HAZARD syntax.
- **[Example vignettes](https://ehrlinger.github.io/temporal_hazard/articles/)** --- exploratory, actuarial, estimation, multivariable, prediction, bootstrap, and sensitivity workflows.

## Development

```r
install.packages(c("devtools", "roxygen2", "pkgdown", "testthat"))
devtools::install_deps(dependencies = TRUE)
devtools::test()
devtools::check()
```

GitHub Actions runs multi-platform `R CMD check` on every push and pull request. Coverage is published to Codecov and the pkgdown site deploys automatically from `main`.
