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
  criterion = c("wald", "aic"), # entry/exit test (see §3)
  slentry = 0.30,               # p-value threshold for entry (criterion = "wald")
  slstay  = 0.20,               # p-value threshold for retention (criterion = "wald")
  max_steps = 50L,
  max_move  = 4L,               # per-variable oscillation limit
  force_in  = character(),      # variables tested but never dropped
  force_out = character(),      # variables never considered
  trace = TRUE,                 # print step-by-step progress
  ...                           # passed to hazard() on refits
)
```

Returns an object of class `c("hzr_stepwise", "hazard")` — the final fit
plus `$steps` (entry/exit trace), `$scope` (candidate set snapshot),
`$criteria` (threshold settings).

**Default thresholds (0.30 entry / 0.20 retention)** match the SAS
`PROC HAZARD` defaults (`SLENTRY = 0.3`, `SLSTAY = 0.2`). The SAS
backward-pass `SLSTAY = 0.05` variant is not modelled separately — one
pair of thresholds governs the two-way loop. Users who want a tighter
backward pass can call `direction = "backward"` with `slstay = 0.05`
after a forward run.

When `criterion = "aic"`, `slentry` / `slstay` are ignored and each step
compares candidate AIC to the current model's AIC: enter if ΔAIC < 0
(lowest ΔAIC wins ties), drop if removing improves AIC.

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
the multiphase scope semantics (per-phase entry, forced-in/out lists)
don't fit the `step()` surface even when `criterion = "aic"`. Keeping a
separate name avoids surprising `stats::step()` users.

---

## 3. Algorithm

### 3.1 Candidate scoring

A single helper `.hzr_candidate_score(fit, variable, phase, criterion)`
returns a numeric scalar with the same "smaller is better" convention
regardless of criterion:

- `criterion = "wald"`: returns the Wald p-value.
- `criterion = "aic"`: returns ΔAIC = AIC_new − AIC_current. Negative
  means the candidate improves the model; positive means it hurts.

Keeping a single scoring abstraction lets the forward / backward /
two-way loops share logic; the only difference is the threshold test at
each step.

#### 3.1a Wald test

For a scalar coefficient β̂_j with SE_j from the fitted vcov:

    z_j = β̂_j / SE_j          p_j = 2·(1 − Φ(|z_j|))

For a multi-df test (e.g. factor with >2 levels), use the joint Wald
statistic:

    W = β̂_S' · V_{S,S}^{−1} · β̂_S  ~  χ²_{|S|}

where S is the index set of the variable's coefficients. This is already
computable from `vcov(fit)`.

- **Retention (backward test):** compute W on the current model — no
  refit.
- **Entry (forward test):** requires a candidate fit with the variable
  added. Exact refit in v1; FAST approximate Wald deferred to v2.

#### 3.1b AIC comparison

    AIC = −2·logLik + 2·k       where k = number of free parameters

For forward entry: fit the candidate, compute
`ΔAIC = AIC_candidate − AIC_current`. Variable enters iff ΔAIC < 0 and
no other candidate has a smaller ΔAIC.

For backward retention: leave-one-out refit is expensive. Use the
single-parameter information-matrix approximation
`ΔAIC ≈ −2·(W_retention) + 2` where W_retention is the Wald χ² statistic
with df = 1 (exact for scalar coefficients; for multi-df variables fall
back to a real refit). Variable drops iff ΔAIC_drop < 0.

### 3.2 Forward step

```
current ← fit
for v in candidates not in current:
  for phase p in v's scope:
    candidate ← hazard(current$formula_extended_by(v, phase), data, dist,
                       phases = current$phases_with(v in phase),
                       theta_start = warm_start(current, v, phase))
    score_v[p] ← .hzr_candidate_score(candidate, v, p, criterion)
(best_v, best_p) ← argmin score_v
threshold_ok ← (criterion = "wald" && score_v[best_v, best_p] < slentry)
             || (criterion = "aic" && score_v[best_v, best_p] < 0)
if threshold_ok:
  add v to phase p
  log step "ENTER v into p  (<criterion-appropriate summary>)"
else: stop
```

The per-candidate refit warm-starts from `current$fit$theta` for the
shared coefficients, with the new coefficient initialised at 0. This
dramatically cuts convergence time for nested fits.

### 3.3 Backward step

```
current ← fit
for v in current model:
  score[v] ← .hzr_candidate_score(current, v, criterion, mode = "drop")
worst ← argmax over v of score[v]  (for Wald — largest p)
         OR argmin over v of score[v]  (for AIC — most negative ΔAIC_drop)
threshold_ok ← (criterion = "wald" && score[worst] > slstay)
             || (criterion = "aic" && score[worst] < 0)
if threshold_ok  AND  worst ∉ force_in:
  drop worst from its phase
  log step "DROP worst from p  (<criterion-appropriate summary>)"
else: stop
```

`force_in` variables are scored and their stats appear in the trace so
the user can see whether they would have dropped, but they are never
actually removed (§2 Q4 decision).

### 3.4 Two-way loop

```
repeat:
  add_happened  ← forward_step()
  drop_happened ← backward_step()
  if !add_happened AND !drop_happened: done
  if max_steps exceeded: warn and stop
  for each variable, track entry+exit count;
    if > max_move: freeze the variable (log "OSCILLATING — frozen")
```

`max_move` prevents a variable that straddles the threshold from cycling
indefinitely. SAS `MOVE` uses the same semantics.

### 3.5 Termination

- `direction = "forward"`: stop when no candidate meets the entry test.
- `direction = "backward"`: stop when no included variable meets the
  drop test (force_in variables are tested but never counted).
- `direction = "both"`: stop when one full pass (forward + backward)
  adds and drops nothing.
- `max_steps` hard cap; emits a warning if hit.
- Non-convergent candidate refits emit a `warning()` naming the failed
  `(variable, phase)` pair and are excluded from the step (§2 Q5
  decision).

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
#                   step_num, action (enter/drop/frozen), variable, phase,
#                   criterion, score, stat, df, logLik, aic, n_coef
#                 `score` is p-value for criterion = "wald", ΔAIC for "aic"
#  $ scope      : list(candidates_per_phase, force_in, force_out)
#  $ criteria   : list(criterion, slentry, slstay,
#                      max_steps, max_move, direction)
#  $ trace_msg  : character vector mirroring console output when trace=TRUE
#  $ final_call : call object for reproducing the final fit
#  $ elapsed    : difftime
#  - attr(*, "class") = c("hzr_stepwise", "hazard")
```

`print.hzr_stepwise()` shows, for `criterion = "wald"`:

```
Stepwise selection (direction = both, criterion = wald,
                    slentry = 0.30, slstay = 0.20)

Step 1: ENTER  nyha    into  early    (p = 0.002)
Step 2: ENTER  age     into  constant (p = 0.018)
Step 3: DROP   shock   from  early    (p = 0.341)
Step 4: ENTER  mal     into  early    (p = 0.103)
(no further action after 4 steps)

Final model:
  <standard hazard summary here>
```

And for `criterion = "aic"`:

```
Stepwise selection (direction = both, criterion = aic)

Step 1: ENTER  nyha    into  early    (ΔAIC = -18.4)
Step 2: ENTER  age     into  constant (ΔAIC =  -6.1)
Step 3: DROP   shock   from  early    (ΔAIC =  -0.9)
(no further action after 3 steps)

Final model:
  <standard hazard summary here>
```

`summary.hzr_stepwise()` is identical to `summary.hazard()` on the final
fit.

---

## 6. Test plan

### 6.1 Unit tests

- `stepwise_direction_forward_wald()`: nested exponential fit where
  adding `x1` is significant; confirm forward selects it.
- `stepwise_direction_forward_aic()`: same scenario under AIC criterion;
  confirm identical variable selected but reported ΔAIC instead of p.
- `stepwise_direction_backward()`: overfitted model; confirm backward
  drops the weakest under each criterion.
- `stepwise_direction_both()`: mixed signal with one oscillator; confirm
  MOVE cap fires and the oscillator is frozen.
- `stepwise_phase_specific()`: multiphase model where `age` is strong in
  `early` but weak in `constant`; confirm it enters only `early`.
- `stepwise_force_in_tested_not_dropped()`: forced variable has
  p = 0.8; confirm it's reported in the trace with its p-value/ΔAIC
  but never dropped (§2 Q4 decision).
- `stepwise_force_out()`: confirm a blocked variable never tries to enter.
- `stepwise_warmstart()`: check that refits converge in ≤ 3 BFGS
  iterations when starting from the prior theta.
- `stepwise_candidate_fails_warns()`: arrange a candidate refit that
  diverges; confirm `expect_warning()` catches a message naming the
  `(variable, phase)` and that the loop continues.
- `stepwise_trace_format()`: snapshot the `steps` tibble structure for
  both criteria.

### 6.2 Parity tests

- **CABGKUL forward selection (Wald)** matches SAS `stepw.c` output
  within numerical tolerance at `SLENTRY = 0.3`. Fixture under
  `inst/extdata/stepwise-fixtures/cabgkul-forward.rds`.
- **AVC multiphase backward (Wald)** matches SAS phase-specific drops
  at `SLSTAY = 0.2`.
- **AIC mode** has no SAS counterpart; validated against `stats::step()`
  on a single-distribution Weibull fit as an external reference.

### 6.3 Edge cases

- All candidates insignificant → forward returns base model unchanged,
  emits a `message()` "no candidates met entry threshold".
- All included variables fail retention → backward returns phase-only
  baseline.
- MOVE = 0 → single-pass forward then single-pass backward (useful for
  deterministic pipelines).
- `direction = "backward"` on a model with no covariates → immediate
  stop with informative message.
- Non-convergence on a candidate refit → emit `warning()` naming
  `(variable, phase)`, skip candidate, continue.

---

## 7. Risks & mitigations

| Risk | Mitigation |
|---|---|
| Refit cost explodes on 30-variable scope × 4 phases | Warm-start from current θ cuts most refits to single BFGS pass; fall back to multi-start only if warm-started fit has non-finite logLik |
| Wald test unreliable for boundary parameters | Reuse existing `summary.hazard()` handling of NA SEs — boundary variables get p = NA and are skipped with a trace note |
| Multiphase phase-specific entry semantics differ from SAS | Write parity tests against `stepw.c` output before claiming parity |
| AIC mode drift from `stats::step()` | Validate single-distribution AIC paths against `stats::step()` output (§6.2) |
| Oscillation despite MOVE cap | Add a hard `max_steps` ceiling and warn on hit |

---

## 8. Implementation order

1. `.hzr_wald_p()` — compute Wald p-value from a fitted `hazard` for a
   named coefficient or coefficient set. Pure function, easy to test.
2. `.hzr_candidate_score()` — unified scoring wrapper dispatching to
   Wald or AIC per the `criterion` argument.
3. `.hzr_phase_update_formula()` / `.hzr_hazard_update_scope()` — formula
   plumbing for add/drop on single-dist and multiphase models.
4. `.hzr_stepwise_forward_step()` — one step of forward, returns updated
   fit or `NULL`.
5. `.hzr_stepwise_backward_step()` — one step of backward.
6. `hzr_stepwise()` — user-facing driver combining the steps with MOVE /
   max_steps / force_in guards.
7. `print.hzr_stepwise()`, `summary.hzr_stepwise()`.
8. Tests 6.1 first, then parity fixtures 6.2.
9. Vignette section "Variable selection" in
   `clinical-analysis-walkthrough.qmd` demonstrating the workflow.

---

## 9. Decisions log

Signed off on 2026-04-16.

| # | Question | Decision |
|---|---|---|
| 1 | API name | `hzr_stepwise()` (no `step.hazard()` S3 method) |
| 2 | Default thresholds | SAS defaults: `slentry = 0.30`, `slstay = 0.20` |
| 3 | AIC mode | Built in from day one via `criterion = c("wald", "aic")` |
| 4 | Force-in semantics | Tested and reported in trace, but never dropped |
| 5 | Non-convergent refit | Emit `warning()` naming `(variable, phase)`, skip, continue |
