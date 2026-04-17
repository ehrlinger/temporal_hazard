# Specify a single hazard phase

Creates an `hzr_phase` object describing one term in a multiphase
additive cumulative hazard model. Pass a list of these to the `phases`
argument of
[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
when `dist = "multiphase"`.

## Usage

``` r
hzr_phase(
  type = c("cdf", "hazard", "constant", "g3"),
  t_half = 1,
  nu = 1,
  m = 0,
  tau = 1,
  gamma = 1,
  alpha = 1,
  eta = 1,
  formula = NULL,
  fixed = character(0)
)

# S3 method for class 'hzr_phase'
print(x, ...)
```

## Arguments

- type:

  Character; phase type: `"cdf"`, `"hazard"`, `"g3"`, or `"constant"`.

- t_half:

  Positive scalar; initial half-life (time at which \\G(t\_{1/2}) =
  0.5\\). Used for `"cdf"` and `"hazard"` phases. SAS early:
  `THALF`/`RHO`.

- nu:

  Numeric scalar; initial time exponent. Used for `"cdf"` and `"hazard"`
  phases. SAS early: `NU`.

- m:

  Numeric scalar; initial shape exponent. Used for `"cdf"` and
  `"hazard"` phases. SAS early: `M`.

- tau:

  Positive scalar; scale parameter for `"g3"` phases. SAS late: `TAU`.

- gamma:

  Positive scalar; time exponent for `"g3"` phases. SAS late: `GAMMA`.

- alpha:

  Non-negative scalar; shape parameter for `"g3"` phases. When
  `alpha > 0`, the generic G3 formula is used; `alpha = 0` gives the
  exponential limiting case. SAS late: `ALPHA`.

- eta:

  Positive scalar; outer exponent for `"g3"` phases. SAS late: `ETA`.

- formula:

  Optional one-sided formula (e.g. `~ age + nyha`) for phase-specific
  covariates. When `NULL` (default), the phase inherits the global
  formula from
  [`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md).

- fixed:

  Character vector naming shape parameters to hold fixed during
  optimization. Valid names for `"cdf"`/`"hazard"`: `"t_half"`, `"nu"`,
  `"m"`, or `"shapes"` (shorthand for all three). Valid names for
  `"g3"`: `"tau"`, `"gamma"`, `"alpha"`, `"eta"`, or `"shapes"`
  (shorthand for all four). Fixed parameters are held at their starting
  values; only `mu` (and covariates) are estimated. Ignored for
  `"constant"` phases. This mirrors the SAS/C HAZARD workflow where
  shapes are typically fixed and only scale parameters are estimated.

- x:

  An `hzr_phase` object (for `print.hzr_phase()`).

- ...:

  Additional arguments (ignored).

## Value

An S3 object of class `"hzr_phase"` with elements:

- type:

  Phase type string.

- t_half:

  Initial half-life (cdf/hazard phases).

- nu:

  Initial time exponent (cdf/hazard phases).

- m:

  Initial shape exponent (cdf/hazard phases).

- tau:

  Scale parameter (g3 phases).

- gamma:

  Time exponent (g3 phases).

- alpha:

  Shape parameter (g3 phases).

- eta:

  Outer exponent (g3 phases).

- formula:

  Phase-specific formula or `NULL`.

- fixed:

  Character vector of fixed parameter names (may be empty).

## Phase types

- `"cdf"`:

  Early risk that resolves over time. \\\Phi(t) = G(t)\\, bounded \\\[0,
  1\]\\. SAS equivalent: Early / G1 phase.

- `"hazard"`:

  Late or aging risk that accumulates. \\\Phi(t) = -\log(1 - G(t))\\,
  monotone increasing. SAS equivalent: Late / G3 phase.

- `"constant"`:

  Flat background hazard rate. \\\Phi(t) = t\\. No shape parameters are
  estimated. SAS equivalent: Constant / G2 phase.

## See also

[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
for fitting multiphase models,
[`hzr_decompos()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_decompos.md)
for the underlying parametric family,
[`hzr_phase_cumhaz()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase_cumhaz.md)
and
[`hzr_phase_hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase_hazard.md)
for computing \\\Phi(t)\\ and \\\phi(t)\\ from these specifications.

[`vignette("fitting-hazard-models")`](https://ehrlinger.github.io/temporal_hazard/articles/fitting-hazard-models.md)
for multiphase fitting examples,
[`vignette("mf-mathematical-foundations")`](https://ehrlinger.github.io/temporal_hazard/articles/mf-mathematical-foundations.md)
for the mathematical framework.

## Examples

``` r
# Classic 3-phase Blackstone pattern
early <- hzr_phase("cdf",      t_half = 0.5, nu = 2, m = 0)
const <- hzr_phase("constant")
late  <- hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1)

# Fix all shapes (C/SAS-style: only estimate mu)
early_fixed <- hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
                          fixed = "shapes")
late_fixed  <- hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1,
                          fixed = "shapes")

# Fix only some parameters
early_partial <- hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
                            fixed = c("nu", "m"))

# Phase with specific covariates
early_cov <- hzr_phase("cdf", t_half = 0.5, nu = 2, m = 0,
                        formula = ~ age + shock)

# Use in hazard():
# hazard(Surv(time, status) ~ age, data = dat,
#        dist = "multiphase",
#        phases = list(early = early, constant = const, late = late))
```
