# TemporalHazard -- Full Package Code Review

**Reviewed:** 2026-04-17
**Package version at review:** 0.9.4 (code changes in response bumped to 0.9.5)
**Scope:** `R/` (21 files, ~9.8 k LOC), `tests/testthat/` (27 files, ~6.8 k LOC),
`man/` (60+ `.Rd` files), 8 Quarto vignettes, `inst/dev/`, `inst/fixtures/`,
`_pkgdown.yml`, `cran-comments.md`. Parity cross-check against
`~/Documents/GitHub/hazard` (the original C/SAS HAZARD codebase).

This document is the canonical write-up of the audit and the follow-up
work it triggered. Cross-reference with:

- `SAS-PARITY-GAP-ANALYSIS.md` -- feature-by-feature parity details
- `DEVELOPMENT-PLAN.md` -- roadmap and phase numbering (Phase 4e/4f
  came directly out of this review)
- `STEPWISE-DESIGN.md` -- v1 scope decisions referenced below
- `NEWS.md` -- version-level scope changes applied in response
- `../../HANDOFF.md` -- session-to-session context transfer

## Overall verdict

Release-candidate quality. Three categories of issue surfaced; all three
share the same pattern -- a feature advertised as complete in NEWS
whose wire-up was incomplete enough to produce silently wrong results
in a narrow regime. Each was resolved by the same recipe: fail fast
with an explicit error, update the NEWS/parity docs to match reality,
add a guard test, track the proper fix in `DEVELOPMENT-PLAN.md`.

Result after the response work (commit `516d12c`):

- `R CMD check --as-cran`: 0 ERROR | 0 WARNING | 2 NOTEs (both cosmetic:
  new submission + HTML Tidy environmental)
- `lintr::lint_package()` clean
- `testthat`: 0 FAIL, 1 WARN, 5 SKIP (all accounted for), 1078 PASS

---

## 1. Documentation completeness

### Fixed during the response

- **Four exported S3 print methods were undocumented.** `print.hazard`,
  `print.summary.hazard`, `print.hzr_phase` were registered in
  `NAMESPACE` via `S3method()` but had no corresponding `.Rd` file. The
  first two conflicted with their parent function's `x` parameter when
  using `@rdname`, so they now carry stand-alone roxygen with
  `@keywords internal` (generates a minimal `.Rd` for R CMD check,
  hidden from pkgdown). `print.hzr_phase` is `@rdname`'d onto
  `hzr_phase.Rd` since `hzr_phase()` takes no `x` argument.
  (`print.hzr_nelson`, `print.hzr_bootstrap`, `print.hzr_competing_risks`
  were already correctly `@rdname`'d onto their parents -- a false
  positive in the original audit.)
- **Stale `cran-comments.md`** (said "TemporalHazard 0.9.1" but DESCRIPTION
  was 0.9.4). Updated to reflect 0.9.5 scope, including the
  `skip_on_cran()` parity tests disclosure.
- **Non-ASCII characters in `R/*.R` triggered a check warning** in
  `R/stepwise-step.R` (em-dash in a user-visible error string on line
  551). Fixed reactively there, then swept systematically across all 21
  R source files: replaced 1429 characters across 33 distinct code
  points (Greek letters -> ASCII names, arrows/math symbols ->
  equivalents, box-drawing characters -> hyphens).
- **Vignette coverage gaps.** `inference-diagnostics.qmd` now demos the
  seven SAS-macro equivalents (`hzr_calibrate`, `hzr_kaplan`,
  `hzr_nelson`, `hzr_gof`, `hzr_deciles`, `hzr_bootstrap`,
  `hzr_competing_risks`) rather than re-implementing each by hand.
  `sas-to-r-migration.qmd` gains a "SAS macro equivalents" section and
  a "Known limitations vs. SAS HAZARD" section naming the five scope
  narrowings honestly.

### Remaining documentation follow-ups

- **Regenerate README figures** once the weighted-CoE narrowing settles.
  `inst/dev/generate-readme-figures.R` exists; re-run before CRAN submit.
- `HANDOFF.md` carries the last session's status; periodic pruning to
  prevent it from becoming a staleness liability.

---

## 2. Feature parity vs. `../hazard` (C/SAS HAZARD)

The original `SAS-PARITY-GAP-ANALYSIS.md` before the review claimed
"Complete" on six of eight capabilities. The audit confirmed five of
them and downgraded two to "Partial". A deep read of the C source also
surfaced **five capabilities that were absent from the gap document
entirely** -- listed at the end of this section.

### Status table (post-review, as of v0.9.5)

| Capability | R status | Notes |
|:---|:---:|:---|
| Multi-phase hazard modeling | Complete | G1/G2/G3 phase types all wired; `hzr_decompos()` implements all six mathematical cases |
| Right + left + interval censoring | Complete | All four `status` codes (-1, 0, 1, 2) handled with correct LL contributions |
| Time-varying covariates | Partial | Global `time_windows` work; SAS allows per-phase EARLY/LATE windows, R applies a single set to all phases (previously unflagged) |
| Phase-specific covariates | Complete | `formula` argument on `hzr_phase()` |
| Weighted events (Weibull, multiphase) | Complete | LL + Weibull analytic gradient correctly weighted (Weibull gradient bug fixed in 0.9.5) |
| Weighted events (exp, log-logistic, log-normal) | **Planned** | Accepted the formal but never applied; `hazard()` now raises an explicit error. Phase 4e. |
| Conservation of Events (unit weights, right-censored) | Complete | Turner's theorem matches SAS to within 1 LL unit on CABGKUL |
| Conservation of Events (non-unit weights) | **Planned** | `.hzr_conserve_events()` sums per-phase cumhaz without weights; auto-disabled when weights != 1 in 0.9.5. Phase 4e. |
| Conservation of Events (interval / left-censored) | **Planned** | Auto-disabled via `coe_supported_data <- all(status %in% c(0, 1))` guard. Phase 4e. |
| Repeating events, trivial case `Surv(0, t, d)` | Complete | Equivalent to `Surv(t, d)`, verified |
| Repeating events, counting-process `Surv(start, stop, event)` with `start > 0` | **Planned** | Parser accepts it but LL never applies `H(stop) - H(start)`. `hazard()` now rejects nonzero starts in 0.9.5. Phase 4f. |
| Stepwise variable selection (forward / backward / two-way) | Complete (v1) | `hzr_stepwise()` with SAS defaults (SLENTRY=0.30, SLSTAY=0.20), per-phase entry, MOVE oscillation guard |
| FAST screening (Lawless-Singhal approx Wald updates) | Planned | v2 enhancement; every candidate currently gets a full refit |
| Multi-df terms in stepwise (factors, splines) | **Planned (v2 by design)** | Upstream guard in `.hzr_candidate_coef_name()` errors with a clear message instructing pre-expansion |
| Vcov / correlation matrix estimation | Complete | `numDeriv::hessian()` at MLE, proper handling of fixed-param NA entries |

### Five gaps the original parity doc did not mention

Surfaced by reading the C source tree:

1. **Prediction confidence limits.** The C tree has a full family
   `hzp_calc_haz_CL.c`, `hzp_calc_srv_CL.c`, `hzp_calc_parm_drv.c`
   that delta-methods CL bands for hazard, survival, and parameter
   derivatives. R `predict.hazard()` returns point estimates only --
   no `se.fit`, `level`, or `*_CL` columns. **Biggest remaining
   functional gap.** Given that `vcov(fit)` is already available, a
   delta-method pass is a reasonable medium-effort addition.
2. **Density / quantile prediction types.** No `type = "density"` or
   `type = "quantile"` (median survival) in `predict.hazard()`.
3. **`OUTEST=` / `OUTVCOV=` dataset export.** Pure ergonomics --
   `coef()` / `vcov()` cover the need in-memory, but SAS veterans
   look for it.
4. **`anova.hazard()` / likelihood-ratio comparison method.** Nothing
   blocks the user from computing -2*dLL by hand, but there's no S3
   method.
5. **`NOCOV` / `NOCOR` summary suppression.** Trivial control flag to
   skip covariance/correlation output; low priority.

All five are documented in `sas-to-r-migration.qmd` -> "Known
limitations vs. SAS HAZARD" so users hit them in the docs rather than
in practice.

---

## 3. Test coverage

`HANDOFF.md` advertised "700+ tests passing" and roughly 62% coverage.
Breadth is good (modules map almost 1:1 to test files) but depth was
uneven.

### Added during the response

- **`test-formula-helpers.R`** (new, 30 tests). Covers right / left /
  interval / counting-process `Surv()` parsing, error paths,
  intercept-only RHS, factor expansion. Previously `R/formula-helpers.R`
  had no dedicated test file.
- **`test-weights.R`** (new, 12 tests). Integer weights equivalent to
  row duplication at LL, gradient, and fit level for Weibull and
  multiphase. Surfaced the Weibull analytic-gradient bug that 0.9.5
  fixes.
- **`test-repeating-events.R`** (new, 7 active + 3 deferred). Parser
  plumbing, trivial-case equivalence, counting-process guard.
  Split-epoch invariance tests `skip("deferred: requires Phase 4f...")`
  so they go live automatically when that feature lands.
- **`test-conservation-of-events.R`** (expanded 7 -> 16). CoE vs
  non-CoE parity, weighted-vs-duplicated fit parity (forces
  `conserve = FALSE`), interval auto-disable, weights auto-disable.

### Modules still without dedicated test files

Low priority (all are covered implicitly through integration tests),
but worth adding for rigour:

- **`R/optimizer.R`** (162 LOC) -- `.hzr_optim_generic()` is the shared
  wrapper behind every distribution; bounds handling, warm-start
  plumbing, and vcov paths are never directly exercised.
- **`R/golden_fixtures.R`** (544 LOC) -- fixture generators are
  assumed-good; no assertion that the `.rds` contents round-trip.
- **`R/argument_mapping.R`** -- `test-argument-mapping.R` exists but is
  26 lines against a 147-line module. Tests the schema shape, not
  the transform rules.

### Parity-test hygiene observations

- `skip_on_cran()` count: 31 in the suite. CRAN sees a lighter test
  set than local dev; disclosed in `cran-comments.md`.
- Parity tolerances in `test-parity-core.R` are generous (1e-3 on
  parameters, 1e-2 on likelihoods). Fine for cross-platform but
  wouldn't catch a subtle CoE regression. Tighten when fixtures are
  regenerated from a stable SAS reference.
- No `testthat` 3rd-edition snapshot tests. Print-method output
  format is regression-exposed. Low priority, but worth adding for
  at least `print.summary.hazard()` and the diagnostic printers.
- `inst/fixtures/hz_loglogistic.rds` and `hz_lognormal.rds` exist on
  disk but are not referenced from any test -- either wire them in or
  delete.

---

## 4. Bugs found and fixed

All three are the same pattern: a feature advertised as complete
whose wire-up missed a code path. Documented so the pattern is
searchable in future audits.

### 4.1 Weibull analytic gradient ignored `weights`

**Symptom.** Numerical gradient on the weighted LL was exactly 2x the
analytic gradient (matching the weighting factor in the test setup).
The optimizer still converged via line-search on the correctly
weighted LL, but the gradient direction was wrong and the final
gradient norm did not approach zero.

**Location.** `R/likelihood-weibull.R`:
- `.hzr_gradient_weibull()` accepted `weights` as a formal but every
  gradient term used `sum(status)` / `sum(cumhaz)` / `crossprod(x, status - H)`
  with no weighting factor.
- `grad_internal()` closure inside `.hzr_optim_weibull()` had the
  same bug plus didn't receive `weights` through its closure at all.

**Fix.** Both gradients now weight the event indicator and cumhaz
building blocks:
```r
w_status <- weights * status
w_cumhaz <- weights * cumhaz
grad[1] <- sum(w_status) / mu - (nu / mu) * sum(w_cumhaz)
# ... etc
```
Verified by `test-weights.R:214` (analytic vs numerical gradient on a
weighted LL, tolerance 1e-5).

### 4.2 Exp / log-logistic / log-normal single-dist likelihoods ignored `weights`

**Symptom.** A weighted fit on any of these distributions silently
returned the unweighted MLE. `weights` was accepted as a formal and
validated, then never referenced inside the LL.

**Location.** `R/likelihood-exponential.R`,
`R/likelihood-loglogistic.R`, `R/likelihood-lognormal.R`. Each has 4
mentions of `weights` -- the formal parameter, the default
assignment, and passthrough in two wrapper calls. Zero mentions in
the actual LL computation.

**Fix (narrow).** `hazard()` now raises an explicit error:
```r
if (!is.null(weights) && !(dist %in% c("weibull", "multiphase"))) {
  stop("'weights' is currently only supported for dist = 'weibull' or ...")
}
```
Full wire-up is tracked in Phase 4e.

### 4.3 Counting-process Surv ignored `time_lower` for status in {0, 1}

**Symptom.** `Surv(start, stop, event)` with `start > 0` scored each
row as `H(stop)` alone instead of `H(stop) - H(start)`. An observation
at `[1.5, 3.0, event=1]` got the same LL contribution as `[0, 3.0,
event=1]`.

**Location.** `.hzr_logl_weibull()` and `.hzr_logl_multiphase()`.
Both compute `cumhaz_lower` but only use it for status == 2
(interval-censored) rows.

**Fix (narrow).** `hazard()` now rejects nonzero-start counting-process
data with an explicit error pointing at the workaround
(`Surv(time, status)` for right-censoring or `Surv(t1, t2,
type="interval2")` for interval). The trivial `Surv(0, t, d)` case
continues to work. Full wire-up is tracked in Phase 4f.

### 4.4 CoE optimizer ignored weights when summing per-phase cumhaz

**Symptom.** Weighted multiphase fit under default settings did not
match the row-duplicated reference fit under the same default
settings, even though the LL path was correct.

**Location.** `.hzr_conserve_events()` in
`R/likelihood-multiphase.R`. Receives the *weighted* event count as
the `total_events` target but then computes `sumcz <-
sum(decomp$total)` and `sumcj <- sum(decomp[[fixmu_phase]])` --
summing raw per-phase cumhaz across rows without weighting. Turner's
adjustment `lfactor <- log(jevent / sumcj)` therefore comes out on a
mismatched scale when weights != 1.

**Fix (narrow).** Added `coe_supported_weights <- all(weights == 1)`
to the CoE guard. Non-unit weights now fall through to the full-
dimensional optimizer, which IS correctly weighted. Fits are still
correct; they just don't benefit from the one-parameter analytical
closed-form solve. Tracked in Phase 4e alongside the other weights
completion tasks.

---

## 5. Side-quest: Copilot PR review

Six Copilot bot comments surfaced on PR #15. Four were stale (anchored
to commit `9eaaee87` and already fixed by `35e8811 "Address Copilot
review on PR #15"` which landed before my audit). Two re-anchored
onto my own commit (`516d12c`) and pointed at the stepwise
main-effects-only v1 scope decision.

Resolved by posting threaded replies on all six:
- 4 stale: "Already fixed at HEAD, comment is anchored to `<old>`"
- 2 live: "Correct observation; v1 scope is deliberately main
  effects only, with an upstream guard that errors clearly on
  multi-df terms. Multi-df support is tracked in
  `STEPWISE-DESIGN.md` as a v2 enhancement"

---

## 6. Remaining action list

Items that survived the narrowing decisions and are genuine future
work. Some are tracked in `DEVELOPMENT-PLAN.md` already; re-listed
here for consolidation.

### Pre-CRAN (blocking)

1. Push CI run for commit `516d12c` should flip `lint` and
   `build-and-deploy` to green (both fixes are in that commit). Verify
   with `gh pr checks 15`.
2. `devtools::spell_check()` -- not yet run.
3. Regenerate README figures via `inst/dev/generate-readme-figures.R`
   after the weighted-CoE narrowing is settled.
4. Run `devtools::check(args = "--as-cran")` one more time on a clean
   build immediately before submission.

### Phase 4e -- Weights wire-up completion (medium effort)

Already broken out in `DEVELOPMENT-PLAN.md` Phase 4e. Summary:

- Multiply `weights[idx_*]` through every LL term in
  `R/likelihood-exponential.R`, `R/likelihood-loglogistic.R`,
  `R/likelihood-lognormal.R`.
- Update each distribution's analytic gradient to use weighted
  building blocks (`w * status`, `w * cumhaz`) -- mirror the 0.9.5
  Weibull fix.
- Thread `weights` from `.hzr_optim_*` into each `grad_internal`
  closure via a captured `.outer_w` (Weibull pattern).
- Remove the `dist %in% c("weibull", "multiphase")` guard in
  `hazard()`.
- Extend `test-weights.R` with per-distribution duplication-parity
  tests.
- Add a `weights` argument to `.hzr_conserve_events()` and apply it
  when summing per-phase cumhaz so the Turner target and
  predicted-event target are on the same scale.
- Remove the `all(weights == 1)` check in the CoE guard
  (`R/likelihood-multiphase.R:951`).

### Phase 4f -- Counting-process LL completion (small-medium effort)

Already broken out in `DEVELOPMENT-PLAN.md` Phase 4f. Summary:

- In `.hzr_logl_weibull()` and `.hzr_logl_multiphase()`, swap
  `cumhaz_event[idx_event] -> cumhaz_event[idx_event] -
  cumhaz_lower[idx_event]` in the event term and the analogous
  subtraction in the right-censored term, honouring `time_lower` for
  all `status` values.
- Analytic gradient: each term picks up an extra `-d/dtheta
  cumhaz_lower` contribution; follows the pattern of the existing
  `cumhaz` derivative code.
- Remove the counting-process-with-start > 0 guard in `hazard()`.
- Remove `skip("deferred: requires Phase 4f ...")` from the three
  disabled tests in `test-repeating-events.R`.
- Add a vignette section on epoch-decomposed longitudinal data.
- Add a SAS parity test against a repeating-events reference fit.

### Unlisted parity gaps

See Section 2 "Five gaps the original parity doc did not mention."
These are not in `DEVELOPMENT-PLAN.md` today.

- **Prediction confidence limits** (delta-method; most valuable).
  Track as a new Phase 4g, or fold into Phase 7 "Performance and
  extensions" if positioned as an enhancement rather than parity.
- Density / quantile prediction types.
- `OUTEST` / `OUTVCOV` dataset export.
- `anova.hazard()` S3 method.
- `NOCOV` / `NOCOR` summary flags.

### Stepwise v2 (scoped for after v1 settles)

Per `STEPWISE-DESIGN.md`:

- Multi-df term support: joint Wald chi-square test,
  column-count-aware `theta` expansion, coefficient-name-matching in
  `.hzr_candidate_coef_name()`.
- FAST screening mode (Lawless-Singhal approximate Wald updates).
- Optional: profile-likelihood confidence intervals for the
  selected-in variables.

### Lower-priority test coverage backfill

- `test-optimizer.R` for `.hzr_optim_generic()` bounds / warm-start /
  vcov paths.
- `test-golden-fixtures.R` for fixture-generator round-trip.
- Expand `test-argument-mapping.R` from 26 lines to cover transform
  rules.
- Add `expect_snapshot()` tests for `print.hazard`, `print.summary.hazard`,
  and the diagnostic print methods.
- Resolve the two orphan fixtures in `inst/fixtures/`
  (`hz_loglogistic.rds`, `hz_lognormal.rds`) -- wire in or delete.
- Tighten parity tolerances when SAS reference fixtures are
  regenerated.

### Spawned separately

- `fix/bootstrap-coef-labels` task: `hzr_bootstrap()$replicates`
  leaves covariate coefficient `parameter` strings blank (only `mu`
  and `nu` get labelled for Weibull). Breaks long-to-wide pivots.
  Scoped as a small standalone fix.

---

## 7. Provenance

This review was conducted on 2026-04-17 across several conversation
phases:

1. **Audit** -- three parallel explore agents against documentation,
   feature parity, and test coverage; synthesis into a four-bucket
   report.
2. **CRAN blockers** -- ASCII sweep, Rd generation for undocumented
   print methods, cran-comments version bump.
3. **Test coverage backfill** (items 3-7 of the audit's "Major" list)
   -- new test files for formula-helpers, weights, repeating events;
   expansion of conservation-of-events; vignette extensions for the
   seven SAS-macro exports + "Known limitations" section.
4. **Bug-finding side effects** -- three narrow-treatment decisions
   for Weibull gradient (fixed), single-dist weights (narrowed),
   counting-process LL (narrowed), CoE + weights (narrowed).
5. **CI follow-through** -- resolved lint + pkgdown failures;
   addressed all six Copilot review comments on PR #15.

Landed as commit `516d12c` on branch `feature/stepwise-selection`.
Sessions logged via the `close` skill in the Obsidian vault under
`Claude/Sessions/2026-04-17 ...`.

---

**Next review trigger.** Phase 4e + 4f completion (the two
narrowings converted to real wire-up) warrants a targeted re-audit
focused on weights and repeating-events coverage, including SAS
parity-test regeneration.
