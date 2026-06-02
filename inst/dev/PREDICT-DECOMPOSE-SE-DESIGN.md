# Per-phase CLs: `predict(decompose=TRUE, se.fit=TRUE)` — Design

**Date:** 2026-06-02
**Roadmap item:** Phase 7c — unblock `predict(decompose = TRUE, se.fit = TRUE)` (DEVELOPMENT-PLAN.md, blocked row)
**Type:** Feature (currently hard-rejected at `hazard_api.R:701`)

## Problem

`predict.hazard()` supports `decompose = TRUE` (per-phase cumulative hazard,
point estimates only) and `se.fit = TRUE` (delta-method CLs on the total)
independently, but the combination is rejected:

    stop("'se.fit = TRUE' with 'decompose = TRUE' is not supported. ...")

There is no per-phase uncertainty quantification. The math to do it already
exists: `.hzr_predict_jacobian_multiphase()` builds each phase's Jacobian block
($\partial H_j/\partial\theta$) before summing them into the total. A per-phase
CL keeps phase $j$'s columns (zeros elsewhere) and runs the same delta-method
sandwich the aggregate path already uses.

## Background (current state)

- `decompose = TRUE, se.fit = FALSE` returns a **wide** data frame
  `data.frame(time, total, <phase1>, <phase2>, ...)`; each phase column is
  $H_j(t) = \mu_j(x)\,\Phi_j(t)$.
- `se.fit = TRUE, decompose = FALSE` returns `data.frame(fit, se.fit, lower,
  upper)` via `.hzr_predict_with_se()`, which: computes `target = diff_fn(theta)`,
  restricts the vcov sandwich to the free-parameter submatrix (handling
  `fixed = "shapes"` / CoE NA rows/cols), subsets the Jacobian columns to match,
  and applies a log-scale back-transform for `cumulative_hazard`.
- Math: $H(t|x) = \sum_j \mu_j(x)\Phi_j(t)$, so each phase's cumulative-hazard
  contribution is additive and its Jacobian block is independent of other
  phases' parameters.

## Decisions

- **Return shape: long.** `data.frame(time, component, fit, se.fit, lower,
  upper)`, one row per time × component, `component ∈ {total, <phase names>}`,
  ordered `c("total", names(phases))`. Tidy / ggplot-native; scales to any phase
  count.
- **Scope: `cumulative_hazard` only.** Per-phase *survival* is non-additive
  (the existing decompose code already punts on per-phase survival), so
  `decompose + se.fit` for `type = "survival"` (and `hazard` /
  `linear_predictor`, already rejected for multiphase) stays a clean error.
- **Per-phase CLs do not sum to the total CL.** Each phase's SE uses only its
  own parameter block; cross-phase covariance contributes only to the total.
  This is correct delta-method behavior and will be documented.
- **No new public arguments.** Reuse `decompose` + `se.fit`.

## Components

### 1. `.hzr_predict_jacobian_multiphase(..., per_phase = FALSE)` — edit (`R/predict-cl.R`)

Today fills phase-$j$ columns into one summed `J` (n × p). Add `per_phase`:

- `per_phase = FALSE` (default): unchanged — returns the summed `n × p` Jacobian.
- `per_phase = TRUE`: returns a **named list** of per-phase `n × p` Jacobians,
  each with only that phase's columns populated (all other columns zero).

DRY: build the per-phase list internally; the summed return becomes
`Reduce(\`+\`, list)`. Column-filling logic is shared between the two modes.

### 2. `.hzr_predict_with_se_decomposed(object, time, x_list, cov_counts, phases, level)` — new (`R/predict-cl.R`)

Mirrors `.hzr_predict_with_se()` but loops over components
`c("total", names(phases))`:

- `total`: summed Jacobian; target $H(t) = \sum_j \mu_j\Phi_j$.
- phase $j$: its per-phase Jacobian; target $H_j(t) = \mu_j\Phi_j$.

Each component reuses the **same** machinery: free-parameter vcov restriction
(`free_idx <- which(is.finite(diag(vcov)))`), `.hzr_predict_se_from_jacobian()`,
and `.hzr_predict_cl_from_se(target, se, level, "log")`. The vcov restriction is
computed once and applied to every component's Jacobian.

Assembles and returns the long data frame; `component` is an ordered factor with
levels `c("total", names(phases))`.

vcov-unavailable / NA-outside-fixed paths return NA CLs (per component) with the
existing warning, matching `.hzr_predict_with_se()` behavior.

### 3. `predict.hazard()` — edit (`R/hazard_api.R`)

Replace the hard `stop()` (currently ~line 701) so that for multiphase models:

- `decompose && se.fit && type == "cumulative_hazard"` → dispatch to component 2
  (wired into the existing multiphase branch ~line 861 where `x_list`/`cov_counts`
  are already assembled).
- `decompose && se.fit && type != "cumulative_hazard"` → clean `stop()`
  (per-phase survival non-additive).

Update the roxygen `@param decompose` / `@param se.fit` / `@details` notes that
currently say the combination is unsupported.

## Data flow

```
predict(obj, newdata, type = "cumulative_hazard",
        decompose = TRUE, se.fit = TRUE, level = 0.95)
  -> multiphase branch assembles pred_time, phases, cov_counts, x_list
  -> .hzr_predict_with_se_decomposed(...)
       per-phase Jacobian list  (.hzr_predict_jacobian_multiphase(per_phase = TRUE))
       total Jacobian = Reduce(`+`, list)
       for component in c("total", phase names):
         J_c[, free_idx]  ->  se = sqrt(diag(J V J^T))  ->  log-scale CL
  -> long data.frame(time, component, fit, se.fit, lower, upper)
```

Example:

```
time  component  fit    se.fit  lower  upper
2.0   total      0.142  0.018   0.110  0.183
2.0   early      0.090  0.012   0.069  0.118
2.0   constant   0.052  0.009   0.037  0.073
5.0   total      ...
```

`nrow = n_time × (1 + n_phases)`; CLs log-scale, so `0 <= lower <= fit <= upper`.

## Error handling

| Condition | Behavior |
|---|---|
| `type = "survival"` (or `hazard`/`linear_predictor`) + `decompose + se.fit` | clean `stop()` (per-phase survival non-additive) |
| vcov unavailable / wrong dim | NA CLs (all components) + existing warning |
| NA in vcov outside fixed rows/cols | NA CLs (all components) + existing warning |
| `fixed = "shapes"` / CoE | per-phase CLs from free-parameter submatrix; no error |
| non-multiphase + `decompose` | unchanged from current behavior |

## Testing (`tests/testthat/test-predict-cl.R`)

1. **Shape:** columns `c(time, component, fit, se.fit, lower, upper)`;
   `levels(component) == c("total", names(phases))`;
   `nrow == n_time * (1 + n_phases)`.
2. **Total regression:** `total` rows equal `predict(decompose = FALSE,
   se.fit = TRUE)` exactly (fit / se.fit / lower / upper).
3. **Per-phase fit:** each phase's `fit` equals `predict(decompose = TRUE,
   se.fit = FALSE)` per-phase cumhaz column.
4. **Analytic vs numeric SE:** per-phase SE matches `numDeriv::jacobian` of an
   $H_j(\theta)$ closure (validates the per-phase Jacobian block).
5. **Ordering invariant:** `lower <= fit <= upper`, all `>= 0`, every component.
6. **Fixed shapes / CoE:** a `fixed = "shapes"` fit yields finite per-phase CLs
   (free-param submatrix path), no error.
7. **Survival still rejected:** `decompose + se.fit, type = "survival"` errors
   cleanly.

## Out of scope

- Per-phase survival CLs (non-additive).
- Wide / list return shapes (long only).
- New public arguments.
- `predict(decompose = TRUE)` for non-multiphase distributions.
