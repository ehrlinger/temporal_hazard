# Hazard Test Fixture Gap List
**Meeting:** John + Rajes  
**Purpose:** Scope the private fixture/sample-job repository that feeds both the SAS/C
and R (TemporalHazard) test harnesses.  
**Context:** v1.0.3 is live on CRAN. This list scopes fixture work for **v1.1+** — no
CRAN urgency applies. All items below are development investments toward a more
complete cross-codebase parity harness.  
**Status of existing work:** 54/54 SAS parity expectations passing against the 8
primary fixtures in `~/Documents/GitHub/hazard/examples/`. The gaps below are
everything not yet covered.

---

## How fixtures are consumed

| Harness | Trigger | Reference file |
|---|---|---|
| R (testthat) | `HAZARD_EXAMPLES_DIR` env var points to examples dir | `test-sas-parity.R`, `test-stepwise-parity.R` |
| C binary | `.sas` + `.lst` paired files | `test-parity-c-binary.R` via `TEMPORAL_HAZARD_BIN` |

New fixtures go in a **private** sibling repo (not `temporal_hazard` or `hazard` public
repos).  R tests skip gracefully when the env vars are unset, so CI on the public
repos is unaffected.

---

## Group A — SAS fixtures exist; R parity tests not yet written

These are quick wins: the reference `.lst` captures are already in `hazard/examples/`,
the parser handles them, and we just need to write the `test_that` blocks.

| Fixture | Feature exercised | Blocker before writing tests |
|---|---|---|
| `hm.death.patient` | Per-patient (not aggregated) data, 12-param multivariable EARLY-phase fit | P1 #1: wide-vcov asymmetry when parameters share names across phases (e.g., `STATUS` appears in both Early and Constant). Fix cursor-based routing first. |
| `hm.death.AVC.deciles` | Multivariable EARLY-phase fit (6 covariates: AGE, COM_IV, MAL, OPMOS, OP_AGE, STATUS) + `%DECILES` → `hzr_deciles()` | Same P1 #1 blocker for vcov; deciles comparison needs `hzr_deciles()` parity tolerance defined |
| `hm.death.AVC` fit 1 | Stepwise selection (`SELECTION SLE=0.2 SLS=0.1`) | Fixture template at `inst/extdata/stepwise-fixtures/cabgkul-forward-wald.sas` exists; must run SAS and call `.hzr_build_stepwise_fixture()` to materialize the `.rds` |
| `hm.dthar.TGA` | Arterial switch (TGA), EARLY-phase only, 9 covariates, E+C, stepwise + final model | `hm.dthar.TGA.sas` created 2026-06-03; needs CCF run to fill `?` PARMS and capture `.lst`. **Public harness** — goes in `hazard/examples/` + `test-sas-parity.R` same as other primary fixtures. |
| `hm.death.AVC` fit 2 | Shaping modifiers `/S`, `/I`, `/E` on covariates | Not in R public API — defer to v1.1 or hand-translate; mark as known gap |
| `ac.death.AVC` | `%KAPLAN` + `%NELSONT` (actuarial Kaplan-Meier and Nelson-Aalen) | Write `hzr_kaplan()` / `hzr_nelson()` parity expectations |
| `hp.death.AVC` | `%HAZPRED` predictions overlaid on KM (hazard, survival, cumulative hazard types) | Define parity tolerance for predict outputs; coordinates with `predict.hazard()` CI tests |
| `hp.death.AVC.hm1`, `hm2` | Pure prediction from saved estimates — two covariate profiles | Same as above |
| `hs.death.AVC.hm1` | Stratified predictions (two strata, two covariate vectors) | `predict.hazard()` newdata path |
| `bs.death.AVC` | `%HAZBOOT(RESAMPL=5)` with embedded stepwise (SLE=0.12, SLS=0.1) | Bootstrap output is non-deterministic; parity metric = variable selection frequency + mean LL, not point estimates |

---

## Group B — R feature exists; need a new SAS/C sample job

These features are tested in R with synthetic data but have no SAS counterpart. We
need to author `.sas` jobs + capture `.lst` output, then wire them into both harnesses.

### B1. Censoring types (beyond standard right-censoring)

| Feature | R test file | SAS HAZARD statement | Notes |
|---|---|---|---|
| Interval censoring (`status = 2`, `Surv(t_low, t_high, event)`) | `test-interval-censoring-weibull.R`, `test-interval-censoring-other-dists.R` | `LCENSOR` + `RCENSOR` pair | Need a dataset where actual intervals exist, not just endpoints. Clinical calibration data (echo re-read windows?) would be ideal. |
| Left censoring (`status = -1`) | `test-formula-helpers.R` (parsing only) | `LCENSOR` only | Verify SAS HAZARD handles pure left-censoring; confirm contribution formula matches. |
| Counting-process / left-truncated at nonzero entry (`Surv(start, stop, event)`) | `test-repeating-events.R` | `LCENSOR STARTTME` (already in `hz.te123.OMC`) | The OMC fixture already covers this! Cross-check the epoch-split invariant test against the OMC expected values to close this. |

### B2. Distribution families (non-multiphase)

| Distribution | R test file | In SAS HAZARD? | Action |
|---|---|---|---|
| Exponential | `test-exponential-dist.R` | Yes (`DIST=EXP` or equivalent) | Create minimal univariate + covariate job |
| Log-logistic | `test-loglogistic-dist.R` | Unknown — verify | If SAS supports it, create sample job. If not, document as R-only. |
| Log-normal | `test-lognormal-dist.R` | Unknown — verify | Same as log-logistic |

### B3. Multiphase structural variants

All existing SAS primary fixtures use either the `cabgkul`, `avc`, `valves`, or `omc`
datasets.  The following phase-structure variants lack dedicated sample jobs:

| Variant | What to exercise | Suggested dataset |
|---|---|---|
| Single-phase Weibull (no multiphase) | Baseline regression; simplest possible SAS HAZARD job | CABGKUL or AVC |
| 2-phase, free shapes on BOTH early AND constant | Currently only AVC has free THALF/NU on early; need both | AVC or TGA |
| 3-phase with at least one free shape | KUL has all-fixed shapes; need free-shape 3-phase | CABGKUL (larger n supports estimation) |
| 3-phase with covariates | No current fixture combines 3-phase + covariate modulation | CABGKUL |
| 2-phase, EARLY-phase covariates only | `hm.death.AVC.deciles` covers this — already in Group A | — |
| 2-phase, BOTH-phase covariates (early + late) | OMC fit 2 has only late-phase; need a job with covariates in both phases | AVC or new synthetic |

### B4. Covariate modulation edge cases

| Edge case | Why it matters | R test status |
|---|---|---|
| Phase-specific covariates, same variable in multiple phases | Tests that `hzr_phase(formula=~x)` applied to both phases produces correct log-likelihood factoring | `test-phase-specific-covariates.R` covers R side; no SAS counterpart |
| Categorical covariate (factor-encoded) | Dummy encoding must match SAS CLASS statement | Not explicitly tested |
| Interaction term in covariate formula | Formula `~ x1 * x2` expansion | Not tested |
| High-dimensional covariate matrix (p > 30) | Optimizer stability, vcov conditioning | `test-parity-edge-cases.R` has p=40 synthetic test; no SAS job |

### B5. Weights

| Scenario | R test | SAS job needed |
|---|---|---|
| Event-only weighting (current OMC morbidity pattern) | `test-weights.R` invariant + `hz.tm123.OMC` | Covered by Group A |
| Uniform non-integer weights (fractional sampling weights) | `test-weights.R` | New job: verify SAS WEIGHT matches R behavior for non-integer weights |
| All-zero-weight rows present (should be excluded) | `test-weights.R` argument validation | Edge case; confirm SAS behavior |

### B6. Conservation of Events (CoE) variants

The P1 #6 gap (Lagrange vs. closed-form CoE for 2-phase early+late) needs a dedicated
fixture to drive the fix:

| Scenario | Status |
|---|---|
| 2-phase E+C (standard): CoE via closed-form μ_C | Covered — `hz.death.AVC` |
| 3-phase E+C+L (standard): CoE distributes across 3 phases | Covered — `hz.deadp.KUL` |
| 2-phase E+L (no constant): CoE needs Lagrange — **P1 #6 gap** | `hz.te123.OMC` fit 1 and `hz.tm123.OMC` document the discrepancy; gap is in the R implementation, not missing fixtures |
| NOCONSERVE explicit: no CoE applied | Covered — `hz.deadp.KUL` (NOCONSERVE) |
| CONSERVE with all phases free-shape | Need a job where CoE interacts with free-shape estimation | No fixture |

---

## Group C — New edge-case fixtures (neither codebase has them)

These need to be authored from scratch: generate or identify appropriate data, write
the SAS job, capture the `.lst`, then write both harness tests.

### C1. Numerical edge cases

| Case | Motivation | Suggested approach |
|---|---|---|
| Very early events (t → 0) | Exercises log-hazard stability near t=0 | Synthetic: Weibull with small scale; real: immediate post-op deaths |
| Very late events (t → ∞) / heavily right-censored data (>90% censored) | Tests tail of survival function, optimizer stability | Synthetic or long-follow-up registry subset |
| Zero events in one phase (phase "empty" of events) | Optimizer must not diverge; μ estimate should hit boundary | Synthetic or carefully selected clinical window |
| Near-duplicate times (many ties) | Affects log-likelihood gradient; potential numerical instability | Synthetic |
| Scale mismatch: event times in days vs. covariates in years | Conditioning of Hessian | Synthetic |

### C2. Sample size extremes

| Case | n | Expected behavior |
|---|---|---|
| Tiny dataset | 20–30 obs, 5–10 events | Optimizer may not converge; vcov singular or ill-conditioned |
| Moderately large | 500–2000 obs (existing datasets are in this range) | Normal operation |
| Large registry | >10,000 obs | Optimizer performance; likelihood evaluation time |

### C3. The TGA dataset

The `tga` dataset is shipped in the package but has no SAS sample jobs. It likely
exercises different mortality patterns than `cabgkul` or `avc`. Worth creating:
- A null-model (no covariates) 2-phase fit
- A multivariable fit using clinically relevant predictors

### C4. Competing endpoints on the same dataset

The OMC data already exercises TE (thromboembolic events) and TM (permanent morbidity)
as paired competing endpoints. The gap is in fully documenting this pattern:

| Pattern | Status |
|---|---|
| Same patient cohort, different event definitions | `hz.te123.OMC` + `hz.tm123.OMC` — covered in Group A |
| Cause-specific hazard (true competing risks, not just two endpoints) | No fixture — out of scope for current SAS HAZARD macro? Confirm with Rajes. |

---

## Group D — Infrastructure decisions for the meeting

These are architectural choices that should be locked down before writing fixtures.

| Decision | Options | Recommendation |
|---|---|---|
| Private fixture repo name and location | (a) `hazard-fixtures` sibling to `hazard` and `temporal_hazard`; (b) subdirectory of `hazard` with `.gitignore` protecting data | Option (a) — cleaner separation, no risk of accidental public exposure |
| Dataset sourcing for new fixtures | (a) existing clinical datasets (CABGKUL, AVC, VALVES, OMC, TGA); (b) synthetic data generated by R | Prefer real data where possible for clinical realism; use synthetic for pure numerical edge cases |
| SAS capture format | `.lst` text files (current pattern) | Retain — parser already handles CRLF + form-feed |
| R fixture format | `.rds` golden files in `inst/fixtures/` (current) vs. raw `.lst` parsed on-the-fly | Retain `.rds` for golden non-SAS fixtures; use live `.lst` parse for SAS parity (current pattern is correct) |
| Env var conventions | `HAZARD_EXAMPLES_DIR` (primary), `HAZARD_OMC_RAW` (OMC raw data), `TEMPORAL_HAZARD_BIN` (C binary) | Document all three in a single `README` in the fixture repo |
| Version tagging | SAS HAZARD version stamped in fixture filenames or README? | README + capture-date comment at top of each `.lst` is sufficient (current pattern) |

---

## Group E — 4-phase models (R-only; no C/SAS parity possible)

The C binary is hard-wired to 3 phases (`phase[1..3]`). TemporalHazard is N-phase
general. These fixtures need R-vs-R round-trip golden files, not `.lst` captures.

The following scenarios are drawn from the vault (`Projects/Temporal Hazard.md`,
`Projects/HAZARD-production-validation.md`):

| Scenario | Phase structure | Dataset | Status | Notes |
|---|---|---|---|---|
| **Acute aortic dissection** | operative + early + constant + late (4 phases) | IRAD published series, CCF production, or Rajes library | ❌ Not started | Primary clinical target for 4-phase. CoE + Case 3 fixes confirmed working on synthetic CABGKUL. First real test of 4-phase in a clinical context. |
| **Hsich / UNOS post-heart-transplant** | 3-phase (early + constant + late) | Eileen Hsich / UNOS registry | ❌ Not started | 3-phase but new dataset. Coordinate with Eileen Hsich; ship de-identified derivative into `data/`. |
| **Sargent CABG reproducibility** | Published model (phase structure from paper) | Sargent CABG paper | ❌ Not started | Extract published LL + MLE table; build a reproducibility vignette. No raw data needed — verify R recovers the published values from the reported starting parameters. |
| **Dual-early (E1 + E2 + C)** | Two CDF early phases + constant | Synthetic or CABGKUL | 🟡 Partially covered | The 4-phase CoE unit tests use this structure synthetically. No vignette or named fixture. |
| **Rajes production library** | 4-phase, phase-specific covariates, shaping modifiers, high-dim, weighted | Blackstone CCF archive | ❌ Externally gated | **Action item: email Rajes** requesting 3–5 production `.sas` + `.lst` pairs. These feed the Production Validation Stream (PV1) and unblock the most clinically realistic fixtures. |

### Sequencing (from vault)

```
Rajes examples → CCF hazard v4.4.6 install → fresh .lst captures → parity fixtures
```

The Rajes email is the single gating dependency for the most important fixtures.
Until those examples arrive, the aortic dissection and Hsich/UNOS fixtures can be
built against synthetic or published data independently.

---

## Prioritization for v1.1+

v1.0.3 is on CRAN. All items below are development investments — priority reflects
feature value and estimated complexity, not release blocking.

| Priority | Item | Estimated effort |
|---|---|---|
| P0 — gates everything else | **Email Rajes** for 3–5 production `.sas`/`.lst` pairs (4-phase, phase-specific covariates, shaping modifiers, high-dim, weighted events) | 15 min |
| P1 — high value, low effort | Run `hm.dthar.TGA.sas` on CCF, fill PARMS, capture `.lst`, add to `hazard/examples/`, write R parity test | 0.5 days (CCF run + test) |
| P1 — high value, low effort | Fix P1 #1 wide-vcov asymmetry; write `hm.death.patient` + `hm.death.AVC.deciles` parity tests | 2–3 days total |
| P1 — high value, low effort | Materialize stepwise parity fixture (`cabgkul-forward-wald.rds`) and confirm `test-stepwise-parity.R` passes | 1 day |
| P1 — high value, low effort | Kaplan + Nelson parity tests (`ac.death.AVC`) | 0.5 days |
| P1 — high value, low effort | Predict parity tests (`hp.*` fixtures) | 1 day |
| P2 — high value, more work | Fix P1 #6 CoE Lagrange for 2-phase E+L; add tighter tolerances to `hz.te123.OMC` fit 1 | 2–3 days |
| P2 — high value, more work | Interval censoring SAS sample job (confirm SAS HAZARD syntax first) | 1–2 days |
| P2 — high value, more work | TGA dataset sample jobs (null model + multivariable) | 1 day |
| P2 — high value, more work | 3-phase free-shape SAS job + R parity | 1 day |
| P2 — high value, more work | 2-phase both-phase covariates SAS job + R parity | 1 day |
| P3 — useful, not urgent | Bootstrap parity (`bs.death.AVC`) | 1 day |
| P3 — useful, not urgent | Log-logistic / lognormal SAS parity (if SAS HAZARD supports) | 0.5 days |
| P3 — useful, not urgent | Numerical edge cases (near-zero times, heavily censored, empty phase) | 1–2 days |
| P3 — useful, not urgent | **Acute aortic dissection 4-phase** R-vs-R fixture + vignette | 1–2 days (data sourcing may gate) |
| P3 — useful, not urgent | **Sargent CABG reproducibility** vignette (published LL + MLEs only) | 1 day |
| P3 — useful, not urgent | **Hsich / UNOS post-HTx** fixture (coordinate with Eileen Hsich) | Externally gated |
| P4 — future / design needed | Shaping modifiers `/S`, `/I`, `/E` in R API (new public API surface) | 3–5 days |
| P4 — future / design needed | Competing-risks / cause-specific hazard | Needs design discussion |

---

## Open questions answered from the hazard codebase

The following were originally listed as questions for Rajes. They are answered definitively
from `src/` — no discussion needed.

---

**Q1: Interval censoring — supported? Syntax?**  
✅ **Yes, fully supported.** Two separate statements:

```sas
LCENSOR starttme;              /* left-truncation / counting-process start */
ICENSOR c3_var = ctime_var;    /* interval censoring: count_var = start_time_var */
```

`ICENSOR` maps to OBS columns 4 (C3 count) and 5 (CTIME). The likelihood in `setlik.c`
uses a Nelson-type approximation: `C3 * ln([CF(T) - CF(CT)] / (T - CT))`. Both
statements are in the yacc grammar (`hazard_y.y`), fully implemented. Parity fixture
is straightforward to write once we have a dataset with actual interval-censored times.

---

**Q2: Log-logistic / lognormal distributions in SAS HAZARD?**  
❌ **No. Not supported.** `structures.h` defines exactly three structs — `early`,
`constant`, `late` — and the entire binary is wired to the three-phase Blackstone/Naftel
G1/constant/G3 model. There is no distribution dispatch, no `DIST=` keyword. These
distributions are R-only extensions in TemporalHazard. Cross-codebase parity is N/A;
those test files stand alone.

---

**Q3: Per-variable modifiers `/S`, `/I`, `/E`, `/O=` (and others) — what do they do?**  
Confirmed from `hazard_y.y` grammar tokens and **EHB institutional knowledge (2026-06-03)**.
These are stepwise selection control modifiers, not covariate transformations.
The `/` syntax uses an undocumented SAS variable-level option feature — **variables
must be comma-separated** when any modifier is used (e.g., `AGE, STATUS/I, COM_IV/S`).

**Per-variable modifiers** (after `/` within EARLY/CONSTANT/LATE blocks):

| Modifier | Token | EHB intent | R equivalent |
|---|---|---|---|
| `/S` | `START` | Warm-start: enter at step 0 but removable. Used when seeding from a prior model and testing interactions — started variables drop out if an added variable subsumes their information. | Not yet in `hzr_stepwise()` |
| `/I` | `INCLUDE` | Force in, never removed. Used for main effects that must stay regardless of p-value (e.g., COMMENCE analysis main effect). | `force_in` ✅ |
| `/E` | `EXCLUDE` | Cannot enter model, but **still generates a 1-step Q-statistic approximation** shown in the "variables not in model" output. Useful for seeing what a variable's effect *would be* without it entering. NOT just documentation — it produces output. | `force_out` 🟡 (R doesn't generate the Q-stat approximation for excluded vars) |
| `/O=N` | `ORDER = NUMBER` | **Essential for interaction analysis.** Controls entry order so main effects always precede their interaction terms. Not cosmetic — drives correct stepwise sequencing. EHB: "essential in interaction analysis." | Not implemented |
| `/MOVE=N` | `MOVE = NUMBER` | Per-variable oscillation cap override | Not implemented |

**SAS syntax note (EHB):** The `/` uses an undocumented SAS variable-level feature.
Comma separation between variables is **required** whenever any modifier is present.
This constraint originates from SAS itself, not from HAZARD.

**Oscillation guard (EHB: "prevent infinite looping when variable goes in→out→in→out"):**
**Already implemented in R** as `max_move` argument to `hzr_stepwise()` (default 4,
in `stepwise.R` lines 7–8). When a variable exceeds `max_move` total moves it is
frozen. The SAS global `MOVE=N` on the `STEPWISE` statement maps to R `max_move`. ✅

**Statement-level stepwise options** (on the `STEPWISE` statement):
`SLENTRY=`, `SLSTAY=`, `MOVE=`, `MAXSTEPS=`, `MAXVARS=`, `STEPWISE` (forward),
`BACKWARD`, `FAST`, `NOPRINTS`, `NOPRINTQ`, `ROBUST`, `SEMIROBUST`, `ONEWAY`.

**Compound variable grouping statements:**

| SAS statement | Meaning | R equivalent |
|---|---|---|
| `RESTRICT varlist` | Forward-only — cannot be removed once entered | Not implemented |
| `BUNDLE <name> varlist` | Variables enter and leave together as a unit | Not implemented |

**Optimizer options in the C binary** (EHB: "3 algorithms"; confirmed from `hazard_y.y` + `setoptim.c`):

| Option | C token | Description |
|---|---|---|
| (default) | — | Quasi-Newton (Dennis & Schnabel); analytical gradient |
| `QUASINEWTON` | `setopt(26)` | Explicit quasi-Newton variant |
| `STEEPEST` | `setopt(27)` | Quasi-Newton with steepest-descent first step — more stable far from optimum |
| `NUMERIC` | `setopt(28)` | Numerical gradients instead of analytical — slower but useful for pathological likelihoods |

**TemporalHazard optimizer** (confirmed from `R/optimizer.R` + `R/likelihood-multiphase.R`):

TemporalHazard uses a three-layer strategy that is substantially more robust than the
C single-start approach:

1. **Multi-start** (`n_starts`, default 5) — random restarts perturbed from user's
   starting θ by `rnorm(sd = 0.5)`; best result kept. C has no multi-start.
2. **Nelder-Mead warmup** (fixed-shape models with 2–10 free params) — runs
   `optim(..., method = "Nelder-Mead", maxit = 5000)` before each BFGS pass.
   Functional analog of C `STEEPEST` — provides a better basin before main optimizer.
3. **Main optimizer:**
   - Unconstrained: `optim(..., method = "BFGS")` with analytical gradient
   - Constrained (Weibull μ, ν > 0): `optim(..., method = "L-BFGS-B")` with lower bounds
4. **Post-fit Hessian** — numerical via `numDeriv::hessian()` by default; analytic via
   `hessian_fn` arg if supplied. C computes analytically during optimization.

**C → R mapping:**

| C option | TemporalHazard equivalent | Notes |
|---|---|---|
| Default quasi-Newton | `BFGS` | Same algorithm family; different implementation |
| `STEEPEST` | Nelder-Mead warmup | Fixed-shape path only |
| `QUASINEWTON` | Always active — no toggle | BFGS is always quasi-Newton |
| `NUMERIC` | Remove `gradient_fn` | Not user-exposed; `optim` falls back to finite-difference |
| Single start | `n_starts = 5` multi-start | Major R advantage: robust to poor starting values |
| `ROBUST` / `SEMIROBUST` | Not yet implemented | Sandwich SEs — Phase 8 candidate |

**Parity implication:** R will match C log-likelihood values when initialized from
the converged C `PARMS` estimates (multi-start then converges to the same basin).
R will often find a *better* optimum from rough starting values because of multi-start.

For the `cabgkul-forward-wald` fixture, the key question for Rajes: does the KUL
example use `/S`, `/O=`, `RESTRICT`, or `BUNDLE`? If only `/I` and plain variables,
`force_in` covers it and the fixture can be materialized now.

---

**Q4: CoE P1 #6 — what does the C binary actually do?**  
The C code (`setcoe_calc_scaling.c`) applies CoE as a **proportional log-scale shift
applied equally to ALL active phase μ parameters**:

```c
lfactr = ln(total_events / sum_cumhaz);   // single scalar
for each active phase j:
    theta[log_mu_j] += lfactr;            // same shift to every phase
```

This is the Lagrange approach: all phases get the same multiplicative factor. R's
current implementation substitutes μ_C (or μ_L in the 2-phase E+L case) in closed
form, leaving the other μ parameters free. That divergence is exactly P1 #6. **The
correct fix is to match the C approach**: compute one scaling factor and apply it
additively to all active log-μ parameters at each optimizer step, rather than
solving for one μ analytically. No new `conserve_method` argument needed —
the C behavior should just become the R behavior.

---

**Q5: TGA dataset — covariates and structure?**  
The R `tga` dataset (470 patients, n=14 variables) matches the SAS flat file exactly.
Available for modeling: `simple`, `dextroin`, `ca_1rl2c`, `hyaaproc`, `no_tca`,
`tca_time`, `age_days`, `arciopyr`, `ca1_2_l`, `opyear`, `source`. Derived
in SAS: `ln_tcat = log(tca_time+1)`, `in_agep1 = 1/((age_days+1)*12/365.2425)`,
`simiagp1 = simple * in_agep1` (interaction).

The existing SAS examples (`hs.dthar.TGA.*`) are **prediction-only** — they consume
a saved model (`TGASW.HMDTHRI1`) that was fit by an earlier job not present in
`hazard/examples/`. A raw HAZARD fit job for TGA **does not exist** in the fixture set
and must be authored. Natural structure: early-phase only (all 9 risk factors from the
SAS comment block), 2-phase E+C or E+C+L depending on the data shape.

---

**Q6: Competing risks / cause-specific hazard in SAS HAZARD?**  
❌ **No native support.** The OBS array has `C1` (events), `C2` (right-censored),
`C3` (interval-censored) — no cause indicator. One event type per PROC HAZARD call.
The OMC TE/TM pattern (two separate PROC HAZARD calls, one per endpoint) **is the
intended design**, not a workaround. Fine-Gray / sub-distribution hazard is out of
scope for both codebases.

---

## Open questions — status after EHB reply (2026-06-03)

**Q1 `/S` modifier scope** — EHB confirmed `/S` is distinct and important (seeding
from prior model + interaction analysis). The fixture question is still whether the
KUL example specifically uses `/S`; ask Rajes. For a first fixture, `/I`-only is
sufficient. `/S` warm-start is a P2 R API addition.

**Q2 TGA HAZARD fit job** — **RESOLVED 2026-06-03.** `tga` flat file confirmed at
`/programs/development/hazard/examples/data/tga` (already in `dist/examples/data/`).
`hm.dthar.TGA.sas` fit job reconstructed and added to `hazard/dist/examples/`.
**TGA goes into the standard public harness** — same as KUL/AVC/VALVES. The `.lst`
output (aggregate statistics, no PHI) will be added to `hazard/examples/` once the
job is run on CCF and the `?` PARMS are filled in. R parity test goes into
`test-sas-parity.R` as a standard Group A test. **Not** a private fixture.

**Q3 Interval censoring** — EHB: used routinely in congenital heart (first operative
day always interval-censored) and in all HCFA/CMS processing. Root cause is fuzzy
patient recall of event dates ("in 1987", "about 5 years ago"). Ask **Carla** whether
she is coding fuzzy dates in the current warehouse. Easy to generate synthetic: set
`ctime = start_date`, `time = start_date + 1 day` (or random small interval). The
C binary handles it natively; the R likelihood is implemented. Blocker is just
finding a real clinical dataset or deciding to go synthetic.

**Remaining open (still need Rajes):**
- Does the KUL forward-wald fixture use `/S`, `/O=`, `RESTRICT`, or `BUNDLE`?
- Where is the TGA flat-file data in the CCF `/programs` archive?
- Sargent CABG paper: where to find it (EHB did not respond to this one)?
