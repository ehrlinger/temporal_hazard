# Hessian Stability at 12+ Parameters — Design

**Status:** Approved design (2026-06-03); pending implementation plan.
**Roadmap row:** DEVELOPMENT-PLAN.md §7c — "Hessian stability at 12+
parameters" (was: *Passing, fragile*).
**Branch model:** all work on `dev` → v1.1.0 (`1.1.0.9000`).

## Problem

Standard errors are computed by inverting a post-fit Hessian. The single
shared chokepoint is `.hzr_optim_generic()` in `R/optimizer.R`, called by
all five distribution families (exponential, Weibull, log-logistic,
log-normal, multiphase). Today it does:

```r
hess <- numDeriv::hessian(objective, par)   # negative log-likelihood
vcov <- solve(hess)                          # bare inversion
```

At 12+ free parameters (multiphase + covariates; e.g. the
`hm.death.AVC.deciles` 13-parameter fit) this is fragile:

- `numDeriv::hessian()` is only symmetric to Richardson tolerance;
  asymmetry alone destabilizes `solve()` at high dimension.
- A near-singular or non-positive-definite Hessian either throws (→ generic
  "not invertible" warning, `NA` SEs) or silently returns garbage / negative
  variances that become `NaN` SEs downstream with no signal.
- No conditioning diagnostic is surfaced, so "passing but borderline" is
  invisible to the user and untested.

## Goals

1. **Robust + loud (universal).** Never silently emit garbage SEs. Detect
   and *name* every degenerate path (ill-conditioned, non-PD, non-finite,
   negative-variance). Use a numerically stable inversion. No user-facing
   API signature change.
2. **More accurate SEs (per family).** Replace the numerical Hessian with an
   analytic Hessian as the production inversion input for all five families,
   validated against numDeriv and against parameterization scale-invariance.

The SAS `.lst` reference for the 13-parameter fit is pursued on a **parallel,
externally-gated track** (PV1 / CCF `hazard` 4.4.6). This design does not
depend on it; validation here is in-package and CRAN-safe.

## Non-goals (adjacent rows, explicitly excluded)

- **P1 #5** — `vcov.hazard()` "invalid 'nrow' value" on fixed-shape fits.
  Separate parser/S3 bug (PRE-CRAN-PARITY-INVENTORY.md).
- **P2 #11** — CoE-projected vcov SE mismatch vs SAS (R returns raw
  inverse-Hessian on the free-parameter submatrix; SAS projects onto the
  conservation manifold). Different row; not conditioning.

These are noted only so the work stays scoped; touching them is out of band.

## Architecture — two layers at the chokepoint

```
.hzr_optim_generic(objective, gradient, theta_start, hessian_fn = NULL, ...)
   │
   ├─ optim() → MLE                              (unchanged)
   ├─ H = hessian_fn(par)        if supplied     ← LAYER 2 (per-family, opt-in)
   │     else numDeriv::hessian(objective, par)  (existing fallback)
   │
   └─ {vcov, rcond, pd} = .hzr_safe_solve(H)     ← LAYER 1 (universal)
```

Both layers operate on the **internal** parameterization, upstream of:
- each family's delta-method `J` transform to natural scale
  (e.g. `R/likelihood-weibull.R` ~L528), and
- the multiphase free→full NA-expansion (`R/likelihood-multiphase.R` ~L1338).

None of that downstream code changes. Approach A (approved): Layer 1 ships
first and universally; Layer 2 rolls out per family, each gated by tests.

## Layer 1 — `.hzr_safe_solve(H)`

New internal helper (in `R/optimizer.R` or a small `R/hessian-invert.R`).
Input: raw Hessian of the negative log-likelihood. Output:
`list(vcov, rcond, pd)`. Steps, in order:

1. **Symmetrize.** `H <- (H + t(H)) / 2`. Always safe; removes
   Richardson-tolerance asymmetry.
2. **Finiteness guard.** If `anyNA(H) || any(!is.finite(H))` →
   `vcov = NA`, warn `"Hessian contains non-finite entries; standard
   errors unavailable"`. Preserves today's NA contract.
3. **Conditioning check.** `rc <- rcond(H)`. If `rc < tol`
   (`tol = .Machine$double.eps^0.5 ≈ 1.5e-8`) → warn, naming the value:
   `"Hessian is ill-conditioned (rcond = <rc>); standard errors may be
   unreliable"`. **Policy: still return the numbers** (do not withhold SEs
   when ill-conditioned) — the user is warned, not denied.
4. **Stable inversion.** Try `chol2inv(chol(H))` (valid + stable iff `H` is
   positive-definite, the expected case at an interior minimum of the NLL;
   `pd <- TRUE`). On `chol()` failure (non-PD) → `pd <- FALSE`, fall back to
   `tryCatch(solve(H), ...)`; if that throws → `vcov = NA` + warning.
5. **Negative-variance guard.** If `any(diag(vcov) < 0)` → warn
   `"Non-positive variance estimates; the optimum may not be a proper
   maximum"`; set offending diagonals and their rows/cols to `NA` (no
   silent `NaN` SEs downstream).
6. **Return** `list(vcov, rcond = rc, pd)`. The chokepoint attaches `rcond`
   and `pd` to the fit result.

**Thresholds (confirmed):** `rcond` warn threshold `1.5e-8`;
warn-but-still-return policy when ill-conditioned.

## Layer 2 — analytic Hessians (production input, per family)

**Plumbing.** `.hzr_optim_generic` gains `hessian_fn = NULL`. When non-NULL,
`H <- hessian_fn(par)` on the internal scale; else numDeriv as today. Each
`.hzr_optim_*` wrapper builds a closure over its data/masks and passes it —
mirroring how `gradient` is already threaded. For multiphase, `hessian_fn`
returns the **free×free** block (respecting `free_mask`, derived on free
params directly to avoid wasted hypergeometric evaluations); the existing
NA-expansion is untouched.

**Coverage contract (fail-loud).** A family's analytic Hessian ships as
production input **only** for the censoring/counting-process/weights
branches where *both* the gradient and the Hessian are analytic. numDeriv
remains the input for any branch not yet converted. Completing remaining
**analytic-gradient** gaps is in-scope as a forcing function; the end-state
target is full analytic coverage, reached incrementally, never half-shipped.

**Rollout order (one PR each, each gated by the §Validation cross-check
before becoming prod input):**

1. **Exponential** — 1 shape param + covariates. Pattern-setter; validates
   the `hessian_fn` plumbing end to end.
2. **Weibull** — 2 params (`alpha`, `psi` internal) + covariates; existing
   delta-method `J` maps to natural scale unchanged.
3. **Log-logistic** — 2 params + covariates; all censoring types +
   counting-process `H(stop) − H(start)` + weights.
4. **Log-normal** — as log-logistic.
5. **Multiphase** — the payload. Second derivatives of the additive
   multiphase NLL w.r.t. `{t_half, ν, m/η, μ}` per phase + covariates,
   including **G3 hypergeometric** second-order shape terms, reusing the
   per-phase machinery the analytic gradient already builds. CoE-constrained
   and fixed-shape fits invert the free-parameter block only.

## Validation (CRAN-safe; no external `.lst`)

Added per family as each analytic Hessian lands:

**(a) Analytic-vs-numDeriv cross-check — the gate.** On a representative
fit, assert `analytic_H ≈ numDeriv::hessian(objective)` elementwise on the
internal scale, relative tolerance `1e-4` (tight Richardson settings on the
numDeriv reference). An analytic Hessian does not become prod input until
this passes. Mirrors the existing analytic-gradient-vs-`numDeriv::grad`
regression pattern.

**(b) Scale-invariance / self-consistency.** Refit with a covariate linearly
rescaled (`x → x/100`), push through the delta-method `J`, assert
natural-scale SEs match the unscaled fit within tolerance. Exercises the
property that is fragile at high dimension.

**(c) 13-parameter anchor.** Dedicated regression on the
`hm.death.AVC.deciles`-style 13-param multiphase fit: converges,
`pd == TRUE`, `rcond` above the warn threshold, all SEs finite and positive,
analytic-vs-numDeriv agree. This flips the §7c row from "passing, fragile"
to "covered."

**(d) Layer-1 unit tests.** Hand-constructed matrices: near-singular
(ill-conditioned warning + still returns), non-PD (Cholesky fallback +
negative-variance NA-guard + warning), non-finite (NA + warning). Pin
Layer 1 independently of any fit.

## Diagnostics surface (additive, minimal)

- `fit$fit$rcond`, `fit$fit$pd` stored from Layer 1.
- `summary.hazard()` prints one line **only when flagged**, e.g.
  `Note: Hessian ill-conditioned (rcond = 3.2e-11); standard errors may be
  unreliable.` Clean fits print nothing new.
- No change to `vcov.hazard()` / `coef()` / `predict()` signatures.

## PR sequencing (all on `dev` → v1.1.0)

1. **PR-1** — Layer 1 + diagnostics + Layer-1 unit tests + 13-param anchor
   (using today's numDeriv Hessian). Closes the actual fragility; may justify
   flipping the row status on its own.
2. **PR-2** — `hessian_fn` plumbing + Exponential analytic Hessian + cross-
   check & invariance tests.
3. **PR-3** — Weibull.
4. **PR-4** — Log-logistic.
5. **PR-5** — Log-normal.
6. **PR-6** — Multiphase analytic Hessian (G3 second-derivative payload) +
   full 13-param cross-check.

Each PR (2–6): close any in-scope analytic-gradient gap, add the analytic
Hessian, gate on the cross-check, add a NEWS entry under the existing
`# TemporalHazard 1.1.0.9000` heading.

## Parity-tolerance handling

When a family's analytic Hessian becomes prod input, re-run the SAS parity
SE asserts (`tolerance = 5e-3`). A shift *toward* the SAS value is expected
(analytic should be more accurate than numDeriv). A shift *away* is a bug
signal. Document per-PR; never silently re-tolerance.
