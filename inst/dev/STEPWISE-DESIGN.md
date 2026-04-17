# Stepwise Covariate Selection — Design Document

**Status:** Design (pre-implementation)
**Target:** Phase 4b of the development plan
**Author:** John Ehrlinger
**Date:** 2026-04-16

This document describes the proposed design for `hzr_stepwise()`, the last
remaining SAS-parity gap. It specifies the user-facing API, algorithm, edge
cases, test plan, and explicit non-goals for the initial implementation.
Sign-off on this design precedes any code.

SAS reference: `src/vars/stepw.c`, `backw.c`, `swvari.c`, `swvarx.c`,
`dfast.c` in the C HAZARD tree.

---

## 1. Scope

### In scope for v1

- Forward, backward, and two-way stepwise selection
- Wald-test-driven entry (`SLENTRY`) and retention (`SLSTAY`) thresholds
- Phase-specific variable entry (a covariate can enter one phase and not
  another)
- MOVE-limited oscillation guard to prevent cycling
- Forced variables (always retained) and blocked variables (never considered)
- Selection trace: step-by-step record of which variables entered / left
  and why
- Tidy output: final `hazard` fit, entry/exit history, and the score table
  evaluated at each step
- Works on all five distributions (`weibull`, `exponential`, `log-logistic`,
  `log-normal`, `multiphase`)

### Out of scope for v1

- **FAST screening** (Lawless–Singhal approximate Wald updates). The full
  refit per candidate is O(p × refits) slower than FAST but is exact and
  sidesteps an entire class of numerical bugs. Add as an optimisation in
  v2 once the reference implementation is validated.
- **Interaction terms** added programmatically during selection (main
  effects only; user can pre-build interactions in their formula).
- **Bootstrap-aggregated selection** — that is `hzr_bootstrap()`'s job and
  already exists.

---

## 2. API

### Primary entry point

```r
hzr_stepwise(
  fit,                          # a fitted `hazard` object (the base model)
  scope = NULL,                 # candidate set (see §2.2)
  direction = c("both", "forward", "backward"),
  slentry = 0.15,               # p-value threshold for entry
  slstay  = 0.15,               # p-value threshold for retention
  max_steps = 50L,
  max_move  = 4L,               # per-variable oscillation limit
  force_in  = character(),      # variables always retained
  force_out = character(),      # variables never considered
  trace = TRUE,                 # print step-by-step progress
  ...                           # passed to hazard() on refits
)
```

Returns an object of class `c("hzr_stepwise", "hazard")` — the final fit
plus `$steps` (entry/exit trace), `$scope` (candidate set snapshot),
`$criteria` (threshold settings).

**Default thresholds (0.15/0.15)** match common clinical-research practice
rather than the SAS defaults (0.3/0.2 forward, 0.05/0.2 backward). SAS
defaults will be offered via `control = list(sas_defaults = TRUE)`.

### Scope specification (§2.2)

Three ways to declare candidates:

1. **Data-driven** (default when `scope = NULL`): all variables in the
   base `data` not already in the model enter the forward-selection pool
   for all phases. Easy and matches `stats::step()` ergonomics.

2. **Per-phase scope** for multiphase models:
   ```r
   scope = list(
     early    = ~ age + nyha + mal + shock,
     constant = ~ age + nyha,
     late     = NULL   # no candidates; phase is baseline-only
   )
   ```
   A variable listed in multiple phases can enter each independently.

3. **Flat scope** for single-distribution models:
   ```r
   scope = ~ age + nyha + mal + shock
   ```

### Method dispatch

Exported as `hzr_stepwise()` only (no `step.hazard()` S3 method). Reason:
`stats::step()`'s contract is based on AIC, not Wald tests, and the
multiphase scope semantics don't fit the `step()` surface. Keeping a
separate name avoids surprising `stats::step()` users.

---

## 3. Algorithm

### 3.1 Wald test for one candidate

For a scalar coefficient β̂_j with SE_j from the fitted vcov:

    z_j = β̂_j / SE_j          p_j = 2·(1 − Φ(|z_j|))

For a multi-df test (e.g. factor with >2 levels), use the joint Wald
statistic:

    W = β̂_S' · V_{S,S}^{−1} · β̂_S  ~  χ²_{|S|}

where S is the index set of the variable's coefficients. This is already
computable from `vcov(fit)` — no new infrastructure.

**Which test for entry vs retention?**

- **Retention (backward test):** trivial — compute W on the fitted model.
- **Entry (forward test):** requires a *candidate* fit with the new
  variable added. This is the expensive step. v1 refits exactly; v2 may
  add FAST approximate Wald.

### 3.2 Forward step

```
current ← fit
for v in candidates not in current:
  for phase p in v's scope:
    candidate ← hazard(current$formula_extended_by(v, phase), data, dist,
                       phases = current$phases_with(v in phase))
    compute p_v at (v, p)
best ← argmin over (v, p) of p_v
if p_best < slentry:
  add v to phase p
  log step "ENTER v into p  (p = ...)"
else: stop
```

The per-candidate refit warm-starts from `current$fit$theta` for the
shared coefficients, with the new coefficient initialised at 0. This
dramatically cuts convergence time for nested fits.

### 3.3 Backward step

```
current ← fit
for v in current model:
  compute retention p_v via Wald test on current vcov
worst ← argmax over v of p_v
if p_worst > slstay  AND  v not in force_in:
  drop v from its phase
  log step "DROP v from p  (p = ...)"
else: stop
```

No refit needed for the *test*; a single refit happens after the drop to
recompute the model without v.

### 3.4 Two-way loop

```
repeat:
  tried_add ← forward_step()
  if tried_add didn't enter: break_add ← TRUE else break_add ← FALSE
  tried_drop ← backward_step()
  if tried_drop didn't drop: break_drop ← TRUE else break_drop ← FALSE
  if break_add AND break_drop: done
  if max_steps exceeded: warn and stop
  for each variable, track entry+exit count;
    if > max_move: freeze the variable (log "OSCILLATING — frozen")
```

`max_move` prevents a variable that straddles the threshold from cycling
indefinitely. SAS `MOVE` uses the same semantics.

### 3.5 Termination

- `direction = "forward"`: stop when no candidate has p < slentry.
- `direction = "backward"`: stop when no included variable has p > slstay.
- `direction = "both"`: stop when one full pass (forward + backward) adds
  and drops nothing.
- `max_steps` hard cap; emits a warning if hit.

---

## 4. Interaction with existing infrastructure

### 4.1 Formula rewriting

For multiphase models, the "formula" is spread across `hzr_phase$formula`
slots. To add variable `v` to phase `early`, rebuild the phase list with
`phases$early$formula <- update(phases$early$formula, ~ . + v)`, then
refit via `hazard()`. A small helper `.hzr_phase_update_formula()` will
hide the plumbing.

### 4.2 Warm-start theta

`hazard()` already accepts `theta_start` (via `control$theta_start` or a
named vector). The stepwise loop passes the current fit's theta expanded
with zeros for new coefficients. Without warm-start, each candidate refit
pays the full multi-start optimization cost (~seconds); with warm-start,
it's typically one BFGS pass (~0.1 s).

### 4.3 CoE interaction

Conservation of Events is a per-fit optimisation. It stays on during
stepwise refits. No special handling required — CoE picks the fixmu phase
per-fit based on current cumhaz contribution.

### 4.4 Interval / left / repeating data

No special path. The Wald test and refit use whatever likelihood the base
fit used. Tested in §6.

---

## 5. Output shape

```r
str(result)
# List of 7
#  $ fit        : hazard object (final model)
#  $ steps      : tibble with columns
#                   step_num, action (enter/drop), variable, phase,
#                   p_value, stat, df, logLik, n_coef
#  $ scope      : list(candidates_per_phase, force_in, force_out)
#  $ criteria   : list(slentry, slstay, max_steps, max_move, direction)
#  $ trace_msg  : character vector mirroring console output when trace=TRUE
#  $ final_call : call object for reproducing the final fit
#  $ elapsed    : difftime
#  - attr(*, "class") = c("hzr_stepwise", "hazard")
```

`print.hzr_stepwise()` shows:

```
Stepwise selection (direction = both, slentry = 0.15, slstay = 0.15)

Step 1: ENTER  nyha    into  early    (p = 0.002)
Step 2: ENTER  age     into  constant (p = 0.018)
Step 3: DROP   shock   from  early    (p = 0.341)
Step 4: ENTER  mal     into  early    (p = 0.103)
(no further action after 4 steps)

Final model:
  <standard hazard summary here>
```

`summary.hzr_stepwise()` is identical to `summary.hazard()` on the final
fit.

---

## 6. Test plan

### 6.1 Unit tests

- `stepwise_direction_forward()`: nested exponential fit where adding `x1`
  is significant; confirm forward selects it.
- `stepwise_direction_backward()`: overfitted model; confirm backward
  drops the weakest.
- `stepwise_direction_both()`: mixed signal with one oscillator; confirm
  MOVE cap fires and the oscillator is frozen.
- `stepwise_phase_specific()`: multiphase model where `age` is strong in
  `early` but weak in `constant`; confirm it enters only `early`.
- `stepwise_force_in()`: confirm a forced variable stays in even with
  p > slstay.
- `stepwise_force_out()`: confirm a blocked variable never tries to enter.
- `stepwise_warmstart()`: check that refits converge in ≤ 3 BFGS
  iterations when starting from the prior theta.
- `stepwise_trace_format()`: snapshot the `steps` tibble structure.

### 6.2 Parity tests

- **CABGKUL forward selection** matches SAS `stepw.c` output within
  numerical tolerance for the same SLENTRY/SLSTAY/dataset. Fixture under
  `inst/extdata/stepwise-fixtures/cabgkul-forward.rds`.
- **AVC multiphase backward** matches SAS phase-specific drops.

### 6.3 Edge cases

- All candidates insignificant → forward returns base model unchanged,
  warning "no candidates met slentry threshold".
- All included variables fail retention → backward returns phase-only
  baseline.
- MOVE = 0 → single-pass forward then single-pass backward (useful for
  deterministic pipelines).
- `direction = "backward"` on a model with no covariates → immediate
  stop with informative message.
- Non-convergence on a candidate refit → skip candidate, log
  `FAILED_REFIT`, continue.

---

## 7. Risks & mitigations

| Risk | Mitigation |
|---|---|
| Refit cost explodes on 30-variable scope × 4 phases | Warm-start from current θ cuts most refits to single BFGS pass; fall back to multi-start only if warm-started fit has non-finite logLik |
| Wald test unreliable for boundary parameters | Reuse existing `summary.hazard()` handling of NA SEs — boundary variables get p = NA and are skipped with a trace note |
| Multiphase phase-specific entry semantics differ from SAS | Write parity tests against `stepw.c` output before claiming parity |
| User expects `step()`-style AIC-based selection | Document prominently; offer `hzr_stepwise_aic()` as v2 if users ask |
| Oscillation despite MOVE cap | Add a hard `max_steps` ceiling and warn on hit |

---

## 8. Implementation order

1. `.hzr_wald_p()` — compute Wald p-value from a fitted `hazard` for a
   named coefficient or coefficient set. Pure function, easy to test.
2. `.hzr_phase_update_formula()` / `.hzr_hazard_update_scope()` — formula
   plumbing for add/drop on single-dist and multiphase models.
3. `.hzr_stepwise_forward_step()` — one step of forward, returns updated
   fit or `NULL`.
4. `.hzr_stepwise_backward_step()` — one step of backward.
5. `hzr_stepwise()` — user-facing driver combining the steps with MOVE /
   max_steps / force_in guards.
6. `print.hzr_stepwise()`, `summary.hzr_stepwise()`.
7. Tests 6.1 first, then parity fixtures 6.2.
8. Vignette section "Variable selection" in
   `clinical-analysis-walkthrough.qmd` demonstrating the workflow.

---

## 9. Open questions for sign-off

1. **API name** — `hzr_stepwise()` as proposed, or align with `stats::step`
   naming via `hzr_step()`?
2. **Default thresholds** — 0.15/0.15 (clinical norm) or SAS's 0.3/0.2?
3. **AIC mode** — skip entirely for v1 (Wald-only), or add a `criterion =
   c("wald", "aic")` switch from the start?
4. **Force-in semantics** — should forced variables also skip the
   retention test entirely (SAS `INCLUDE`) or be tested but never dropped
   regardless of p?
5. **Boundary / non-convergent candidates** — skip silently with a trace
   note (current proposal) or surface as a warning?
