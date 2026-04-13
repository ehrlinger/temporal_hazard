# Package index

## Hazard model

- [`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
  : Build and optionally fit a hazard model
- [`summary(`*`<hazard>`*`)`](https://ehrlinger.github.io/temporal_hazard/reference/summary.hazard.md)
  : Summarize a hazard model
- [`predict(`*`<hazard>`*`)`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
  : Predict from a hazard model object
- [`coef(`*`<hazard>`*`)`](https://ehrlinger.github.io/temporal_hazard/reference/coef.hazard.md)
  : Extract coefficients from hazard model
- [`vcov(`*`<hazard>`*`)`](https://ehrlinger.github.io/temporal_hazard/reference/vcov.hazard.md)
  : Extract variance-covariance matrix from hazard model

## Multiphase phases

- [`hzr_phase()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase.md)
  : Specify a single hazard phase
- [`is_hzr_phase()`](https://ehrlinger.github.io/temporal_hazard/reference/is_hzr_phase.md)
  : Test if an object is an hzr_phase
- [`hzr_phase_cumhaz()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase_cumhaz.md)
  : Cumulative hazard contribution from a single phase
- [`hzr_phase_hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase_hazard.md)
  : Instantaneous hazard contribution from a single phase

## Decomposition

- [`hzr_decompos()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_decompos.md)
  : Generalized temporal decomposition
- [`hzr_decompos_g3()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_decompos_g3.md)
  : Late-phase (G3) temporal decomposition

## Migration tools

- [`hzr_argument_mapping()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_argument_mapping.md)
  : Legacy HAZARD to hvtiRhazard argument mapping

## Numerical Primitives

- [`hzr_log1pexp()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_log1pexp.md)
  : Numerically stable log(1 + exp(x))
- [`hzr_log1mexp()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_log1mexp.md)
  : Numerically stable log(1 - exp(-x)) for x \> 0
- [`hzr_clamp_prob()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_clamp_prob.md)
  : Clamp probabilities away from 0 and 1
