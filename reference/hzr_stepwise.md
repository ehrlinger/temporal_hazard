# Stepwise covariate selection for a parametric hazard model

Run forward, backward, or two-way stepwise selection on an existing
`hazard` fit using Wald p-values or AIC deltas as the entry / retention
criterion. Phase-specific entry is supported for multiphase models: a
covariate can enter one phase and not another.

## Usage

``` r
hzr_stepwise(
  fit,
  scope = NULL,
  data,
  direction = c("both", "forward", "backward"),
  criterion = c("wald", "aic"),
  slentry = 0.3,
  slstay = 0.2,
  max_steps = 50L,
  max_move = 4L,
  force_in = character(),
  force_out = character(),
  trace = TRUE,
  ...
)
```

## Arguments

- fit:

  A fitted `hazard` object built via the
  `formula = Surv(...) ~ predictors, data = df` interface.

- scope:

  Candidate set. `NULL` (default) uses every data-frame column not
  already in the model for every phase. For single-distribution fits,
  pass a one-sided formula (`~ age + nyha`) or a character vector of
  names. For multiphase fits, pass a named list of one-sided formulas
  keyed by phase.

- data:

  Data frame the base fit was built on. Required for refits.

- direction:

  One of `"both"` (default), `"forward"`, `"backward"`.

- criterion:

  One of `"wald"` (default) or `"aic"`. SAS-style p-value thresholds
  apply to Wald; AIC uses `DeltaAIC < 0` uniformly.

- slentry:

  Entry p-value threshold for the Wald criterion. Default `0.30` matches
  SAS `SLENTRY`.

- slstay:

  Retention p-value threshold for the Wald criterion. Default `0.20`
  matches SAS `SLSTAY`.

- max_steps:

  Hard cap on total accepted actions. Emits a
  [`warning()`](https://rdrr.io/r/base/warning.html) if hit. Default
  `50`.

- max_move:

  Per-variable oscillation cap. When a variable has entered + exited
  more than `max_move` times it is frozen for the remainder of the run.
  Default `4`.

- force_in:

  Character vector of variables that must remain in the model. Such
  variables are still scored and reported in the selection trace, but
  are never dropped.

- force_out:

  Character vector of variables that may never be considered as
  candidates.

- trace:

  Logical; print step-by-step progress to the console. Default `TRUE`.

- ...:

  Passed to the underlying
  [`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
  refits (e.g. `control = list(n_starts = 3)`).

## Value

An object of class `c("hzr_stepwise", "hazard")` – the final fit
augmented with:

- `steps`:

  Data frame with one row per accepted / frozen action; see Details.

- `scope`:

  Record of the candidate scope, plus `force_in`, `force_out`, and the
  frozen set.

- `criteria`:

  Named list of the threshold / direction settings actually applied.

- `trace_msg`:

  Character vector of the trace lines, captured regardless of the
  `trace` flag.

- `elapsed`:

  `difftime` from start to finish.

- `final_call`:

  The call that produced this result.

## Details

The `steps` data frame has columns:

- `step_num`:

  Integer sequence starting at 1.

- `action`:

  `"enter"`, `"drop"`, or `"frozen"`.

- `variable`:

  Variable affected.

- `phase`:

  Phase name (multiphase) or `NA_character_`.

- `criterion`:

  `"wald"` or `"aic"`.

- `score`:

  Winning score used for the decision.

- `stat`, `df`:

  Wald statistic and degrees of freedom.

- `p_value`, `delta_aic`:

  Always populated when computable, regardless of the active criterion.

- `logLik`, `aic`, `n_coef`:

  Goodness-of-fit diagnostics of the model *after* this step.
