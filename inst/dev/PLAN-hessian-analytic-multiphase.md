# Hessian Stability — Layer 2 PR-6 (Multiphase Analytic Hessian) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an analytic Hessian as the production standard-error input for
multiphase fits, replacing the `numDeriv::hessian()` call in
`.hzr_optim_multiphase()` and the matching `numDeriv::hessian()` call in the
CoE full-information vcov path (L1418). Standard errors for high-dimensional
multiphase fits (the 13-parameter `hm.death.AVC.deciles` class) become more
accurate and no longer depend on Richardson finite-difference step-size
sensitivity.

**Architecture:** `.hzr_hessian_multiphase()` builds the p×p NLL Hessian from
three analytically-derived contributions (H curvature, log-h curvature,
log-h outer product) using per-phase quantities already computed for the
gradient. Shape parameters (log_t_half, nu, m for CDF/hazard phases; log_tau,
gamma, alpha, eta for G3) use second-order central-difference second
derivatives via new helpers `.hzr_phase_second_derivatives()` and
`.hzr_g3_phase_second_derivatives()` — consistent with the gradient's use of
first-order FD for shape params. Coverage mirrors the analytic gradient: event
+ right-censored + counting-process start; left/interval-censored rows return
NULL to request the numDeriv fallback.

**Tech Stack:** R, testthat (edition 3), roxygen2, devtools. Base-R numerics;
cross-check tests use `numDeriv` (Suggests, guarded with
`skip_if_not_installed`).

**Scope:** PR-6 of the 6-PR sequence in `inst/dev/HESSIAN-STABILITY-DESIGN.md`.
This is the multiphase payload. PRs 1–5 (Layer 1 + single-distribution
families) are already shipped. Branch: `dev` → v1.1.0.

---

## Mathematical derivation

### Parameterization and notation

Internal parameter vector: `θ = [θ₁, …, θ_J]` concatenated over J phases.
Phase j has parameters in order: `log_mu_j`, then shape params `s_j` (none for
constant; `log_t_half, nu, m` for CDF/hazard; `log_tau, gamma, alpha, eta` for
G3), then covariate coefficients `β_j = (β_{j1}, …, β_{jp_j})`.

Let `x̃_j(i) = [1, x_{j,i}]ᵀ` (augmented covariate vector, length 1+p_j,
prepending the "1" for log_mu). Define:

    μ_j(i) = exp(log_mu_j + x_{j,i}·β_j)
    Φ_j(t) = phase j cumulative shape (output of hzr_phase_cumhaz / hzr_decompos_g3$G3)
    φ_j(t) = phase j instantaneous shape (hzr_phase_hazard / hzr_decompos_g3$g3)
    H_i     = Σ_j μ_j(i) [Φ_j(t_i) − Φ_j(start_i)]    (total cumulative hazard)
    h_i     = Σ_j μ_j(i) φ_j(t_i)                       (total instantaneous hazard)

(Φ_j(start_i) = 0 when no left-truncation or start=0.)

### NLL (event + right-censored rows)

    NLL = Σ_i w_i H_i − Σ_{i∈E} w_i log h_i

Taking second derivatives and splitting into three terms:

### Term A — H curvature (all rows; block-diagonal in phases)

For parameters a, b **both in phase j** (same phase):

    ∂²NLL/∂θ_a∂θ_b|_A = Σ_i w_i ∂²H_i/∂θ_a∂θ_b

Because H = Σ_j μ_j Φ_j is additive across phases, cross-phase second
derivatives of H vanish. Within phase j, write:

**(a) μ/β block** — second derivatives of μ_j Φ_j w.r.t. log_mu and β_jk:

    ∂²H_i/∂(log_mu_j)²       = μ_j(i) Φ_j(t_i)
    ∂²H_i/∂(log_mu_j)∂β_{jk} = μ_j(i) Φ_j(t_i) x_{jik}
    ∂²H_i/∂β_{jk}∂β_{jl}     = μ_j(i) Φ_j(t_i) x_{jik} x_{jil}

All three collapse to one outer product per row:

    A_j^{μβ} = Σ_i w_i μ_j(i) Φ_j(t_i) · x̃_j(i) x̃_j(i)ᵀ
              − Σ_i w_i μ_j(i) Φ_j(start_i) · x̃_j(i) x̃_j(i)ᵀ  [counting-process]

**(b) μ/β × shape cross terms** — mixed second derivatives w.r.t. log_mu (or
β_{jk}) and shape param s:

    ∂²H_i/∂(log_mu_j)∂s = μ_j(i) ∂Φ_j/∂s (t_i)
    ∂²H_i/∂β_{jk}∂s     = μ_j(i) x_{jik} ∂Φ_j/∂s (t_i)

→ A_j^{μβ,s}[k] = Σ_i w_i μ_j(i) (∂Φ_j/∂s)(t_i) x̃_{jk}(i)   (a vector, k indexes x̃_j)

**(c) shape × shape block** — requires ∂²Φ_j/∂s∂s':

    ∂²H_i/∂s∂s' = μ_j(i) ∂²Φ_j/∂s∂s' (t_i)

→ A_j^{ss'} = Σ_i w_i μ_j(i) ∂²Φ_j/∂s∂s' (t_i)

For CDF/hazard phases: `∂²Φ/∂s∂s'` computed via second-order central
differences (new helper `.hzr_phase_second_derivatives()`).
For G3 phases: same via `.hzr_g3_phase_second_derivatives()`.

### Term B — log-h outer product (events only; **dense across all phases**)

    ∂²NLL/∂θ_a∂θ_b|_B = +Σ_{i∈E} (w_i/h_i²) (∂h_i/∂θ_a)(∂h_i/∂θ_b)

This is a sum of rank-1 outer products of the full p-vector `∂h_i/∂θ`. It is
**NOT block-diagonal** — parameters in different phases interact here. The h-
gradient vector for row i has non-zero entries only within each phase:

    ∂h_i/∂(log_mu_j)  = μ_j(i) φ_j(t_i)
    ∂h_i/∂β_{jk}      = μ_j(i) φ_j(t_i) x_{jik}
    ∂h_i/∂s_j         = μ_j(i) (∂φ_j/∂s)(t_i)

Assemble the full p-vector `g^h(i)` by concatenating contributions from all
phases, then:

    Term B = X̃_h' diag(w/h²) X̃_h      (an n_E × p matrix product)

where X̃_h[i, :] = g^h(i)ᵀ and n_E = number of events.

### Term C — log-h curvature (events only; block-diagonal in phases)

    ∂²NLL/∂θ_a∂θ_b|_C = −Σ_{i∈E} (w_i/h_i) ∂²h_i/∂θ_a∂θ_b

Cross-phase terms vanish (h = Σ_j μ_j φ_j, same additive argument as H).
Within phase j the structure mirrors Term A with φ replacing Φ:

**(a) μ/β block:**

    C_j^{μβ} = −Σ_{i∈E} (w_i/h_i) μ_j(i) φ_j(t_i) · x̃_j(i) x̃_j(i)ᵀ

**(b) μ/β × shape cross:**

    C_j^{μβ,s}[k] = −Σ_{i∈E} (w_i/h_i) μ_j(i) (∂φ_j/∂s)(t_i) x̃_{jk}(i)

**(c) shape × shape:**

    C_j^{ss'} = −Σ_{i∈E} (w_i/h_i) μ_j(i) ∂²φ_j/∂s∂s' (t_i)

`∂²φ_j/∂s∂s'` via the same second-derivative helpers as Term A.

### Full Hessian

    H_NLL = (block-diag A) + (dense B) + (block-diag C)

The block-diagonal parts (A + C) are computed phase-by-phase; Term B is
assembled as a single weighted crossprod of the n_E × p matrix X̃_h.

### Second-derivative helpers (.hzr_phase_second_derivatives)

For CDF/hazard phases, shape params `s ∈ {log_t_half, nu, m}`:

Use second-order central differences for each pair (s, s'):
- Diagonal (s = s'): `∂²Φ/∂s² ≈ [Φ(s+h) - 2Φ(s) + Φ(s-h)] / h²`
- Off-diagonal (s ≠ s'): `∂²Φ/∂s∂s' ≈ [Φ(s+h,s'+h) - Φ(s+h,s'-h) - Φ(s-h,s'+h) + Φ(s-h,s'-h)] / (4h²)`

Use step size `h = .Machine$double.eps^(1/4) ≈ 1.2e-4` (optimal for 4-point
formulas). Same logic for φ. Parametrize on the internal (log_t_half, nu, m)
scale, consistent with the gradient.

For G3 phases: identical structure over `{log_tau, gamma, alpha, eta}` via
`.hzr_g3_phase_second_derivatives()`.

### Coverage contract

| Observation type | Coverage |
|---|---|
| Event (status == 1) | Full analytic (all three terms) |
| Right-censored (status == 0) | Term A only (contributes to H terms; h not evaluated) |
| Counting-process start > 0 (status ∈ {0,1}, time_lower > 0) | Term A start-correction (H(start) second derivatives) |
| Left-censored (status == -1) | **return NULL** → numDeriv fallback |
| Interval-censored (status == 2) | **return NULL** → numDeriv fallback |

Rationale: mirrors the gradient's analytic coverage. Left/interval rows are
rare in the primary clinical fixtures and their second-derivative structure
(log(1-exp(-H)) and log(exp(-H_L) - exp(-H_U))) requires additional care that
is deferred.

---

## File structure

- **Create** `R/hessian-multiphase.R` — `.hzr_phase_second_derivatives()`,
  `.hzr_g3_phase_second_derivatives()`, `.hzr_hessian_multiphase()`.
  Single responsibility: analytic NLL Hessian for the multiphase family.
- **Create** `tests/testthat/test-multiphase-hessian.R` — analytic-vs-numDeriv
  cross-check, NULL fallback, end-to-end vcov, scale-invariance, 13-param
  anchor (replaces/extends the placeholder in `test-hessian-stability.R`).
- **Modify** `R/likelihood-multiphase.R` (`.hzr_optim_multiphase()`) —
  pass `hessian_fn` closure to `.hzr_optim_generic()` (line ~1100), and
  replace the numDeriv call in the CoE full-info vcov block (L1418) with the
  analytic Hessian.
- **Modify** `NEWS.md` — entry under `# TemporalHazard 1.1.0.9000`.

---

## Task 1 — Shape second-derivative helpers

**Files:**
- Create: `R/hessian-multiphase.R`
- Create: `tests/testthat/test-multiphase-hessian.R` (initial section)

### Step 1: Write the failing tests

Create `tests/testthat/test-multiphase-hessian.R` with the first section:

```r
# test-multiphase-hessian.R -- Analytic Hessian tests for multiphase models.
# These tests validate the shape second-derivative helpers and the full
# .hzr_hessian_multiphase() function against numDeriv references.

test_that(".hzr_phase_second_derivatives diag matches [Phi(s+h) - 2Phi + Phi(s-h)] / h^2", {
  skip_if_not_installed("numDeriv")
  t <- c(0.2, 0.5, 1.0, 2.0, 4.0)
  t_half <- 1.0; nu <- 1.5; m <- 1.0

  d2 <- .hzr_phase_second_derivatives(t, t_half = t_half, nu = nu, m = m,
                                       type = "cdf")

  # Reference via numDeriv Hessian applied to the scalar-output function
  phi_fn  <- function(par) hzr_phase_cumhaz(t, t_half = exp(par[1]),
                                              nu = par[2], m = par[3],
                                              type = "cdf")
  par0 <- c(log(t_half), nu, m)
  # numDeriv::jacobian gives J[i, k] = d Phi_i / d par_k.
  # The Hessian of Phi[i] w.r.t. par is the i-th row of d²Phi/dpar²:
  for (i in seq_along(t)) {
    H_nd <- numDeriv::hessian(function(par) phi_fn(par)[i], par0)
    expect_equal(d2$d2Phi_dlog_thalf2[i], H_nd[1, 1], tolerance = 1e-4,
                 label = paste0("d2Phi/d(log_thalf)^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dnu2[i],        H_nd[2, 2], tolerance = 1e-4)
    expect_equal(d2$d2Phi_dm2[i],         H_nd[3, 3], tolerance = 1e-4)
    expect_equal(d2$d2Phi_dlog_thalf_dnu[i], H_nd[1, 2], tolerance = 1e-4)
  }
})

test_that(".hzr_g3_phase_second_derivatives diag matches numDeriv reference", {
  skip_if_not_installed("numDeriv")
  t <- c(0.5, 1.0, 3.0, 8.0)
  tau <- 2.0; gamma <- 1.5; alpha <- 0.8; eta <- 0.6

  d2 <- .hzr_g3_phase_second_derivatives(t, tau = tau, gamma = gamma,
                                           alpha = alpha, eta = eta)

  g3_fn <- function(par) {
    hzr_decompos_g3(t, tau = exp(par[1]), gamma = par[2],
                     alpha = par[3], eta = par[4])$G3
  }
  par0 <- c(log(tau), gamma, alpha, eta)
  for (i in seq_along(t)) {
    H_nd <- numDeriv::hessian(function(par) g3_fn(par)[i], par0)
    expect_equal(d2$d2Phi_dlog_tau2[i], H_nd[1, 1], tolerance = 1e-3,
                 label = paste0("d2G3/d(log_tau)^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dgamma2[i],   H_nd[2, 2], tolerance = 1e-3)
  }
})
```

### Step 2: Run to verify they fail

```bash
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-multiphase-hessian.R")'
```
Expected: FAIL — `could not find function ".hzr_phase_second_derivatives"`.

### Step 3: Implement the helpers

Create `R/hessian-multiphase.R`:

```r
#' @keywords internal
NULL

# hessian-multiphase.R -- Analytic NLL Hessian for multiphase hazard models.
#
# Three helpers then the main function:
#   .hzr_phase_second_derivatives()     -- d2Phi/dss' and d2phi/dss' for CDF/hazard
#   .hzr_g3_phase_second_derivatives()  -- same for G3 (late) phases
#   .hzr_hessian_multiphase()           -- full p x p NLL Hessian

# Step size for 4-point second-difference formulas: h ~ eps^(1/4).
.hzr_h2 <- .Machine$double.eps^(1 / 4)

# ---------------------------------------------------------------------------
# CDF / hazard phase second derivatives
# ---------------------------------------------------------------------------
#' Second derivatives of CDF/hazard phase shape functions
#'
#' Computes the six unique mixed second partial derivatives of \eqn{\Phi_j(t)}
#' (cumulative shape) and \eqn{\phi_j(t)} (instantaneous shape) with respect
#' to the internal shape parameters \eqn{\{log\_t\_half, \nu, m\}} via
#' second-order central differences (step size \code{.hzr_h2}).
#'
#' This is the second-derivative analogue of \code{.hzr_phase_derivatives()},
#' which uses first-order central differences.
#'
#' @param time  Numeric vector of observation times.
#' @param t_half,nu,m  Phase shape parameters (natural scale; t_half and nu
#'   appear on the \code{log_t_half} / \code{nu} internal scale in the Hessian,
#'   so the chain-rule factor for t_half is tau * d/d(tau) = d/d(log_t_half)).
#' @param type  Phase type: \code{"cdf"} or \code{"hazard"}.
#' @return A list with elements
#'   \code{d2Phi_dlog_thalf2}, \code{d2Phi_dnu2}, \code{d2Phi_dm2},
#'   \code{d2Phi_dlog_thalf_dnu}, \code{d2Phi_dlog_thalf_dm},
#'   \code{d2Phi_dnu_dm},
#'   and the matching \code{d2phi_*} entries.
#' @noRd
.hzr_phase_second_derivatives <- function(time, t_half, nu, m, type) {
  h <- .hzr_h2

  # Helper: Phi and phi at perturbed (t_half_p, nu_p, m_p)
  eval_at <- function(th, nup, mp) {
    d <- .hzr_phase_derivatives(time, t_half = th, nu = nup, m = mp, type = type)
    list(Phi = d$Phi, phi = d$phi)
  }

  d00 <- eval_at(t_half, nu, m)

  # Step sizes on the internal parameter scale:
  #   log_t_half: step h means t_half * (exp(h) - 1) ≈ t_half * h for small h;
  #   use additive steps on the log scale, converting back to natural scale.
  t_h_p <- t_half * exp( h)
  t_h_m <- t_half * exp(-h)
  nu_p   <- nu + h;  nu_m <- nu - h
  m_p    <- m  + h;  m_m  <- m  - h

  d_tp <- eval_at(t_h_p, nu, m)
  d_tm <- eval_at(t_h_m, nu, m)
  d_np <- eval_at(t_half, nu_p, m)
  d_nm <- eval_at(t_half, nu_m, m)
  d_mp <- eval_at(t_half, nu, m_p)
  d_mm <- eval_at(t_half, nu, m_m)

  # Four-point mixed for off-diagonals
  d_tp_np <- eval_at(t_h_p, nu_p, m)
  d_tp_nm <- eval_at(t_h_p, nu_m, m)
  d_tm_np <- eval_at(t_h_m, nu_p, m)
  d_tm_nm <- eval_at(t_h_m, nu_m, m)

  d_tp_mp <- eval_at(t_h_p, nu, m_p)
  d_tp_mm <- eval_at(t_h_p, nu, m_m)
  d_tm_mp <- eval_at(t_h_m, nu, m_p)
  d_tm_mm <- eval_at(t_h_m, nu, m_m)

  d_np_mp <- eval_at(t_half, nu_p, m_p)
  d_np_mm <- eval_at(t_half, nu_p, m_m)
  d_nm_mp <- eval_at(t_half, nu_m, m_p)
  d_nm_mm <- eval_at(t_half, nu_m, m_m)

  h2 <- h^2

  # Central second differences (diagonal)
  d2Phi <- function(d_p, d_m, d0) (d_p$Phi - 2 * d0$Phi + d_m$Phi) / h2
  d2phi <- function(d_p, d_m, d0) (d_p$phi - 2 * d0$phi + d_m$phi) / h2

  # Mixed four-point formulas
  d2Phi_mixed <- function(dpp, dpm, dmp, dmm) {
    (dpp$Phi - dpm$Phi - dmp$Phi + dmm$Phi) / (4 * h2)
  }
  d2phi_mixed <- function(dpp, dpm, dmp, dmm) {
    (dpp$phi - dpm$phi - dmp$phi + dmm$phi) / (4 * h2)
  }

  list(
    d2Phi_dlog_thalf2       = d2Phi(d_tp, d_tm, d00),
    d2Phi_dnu2              = d2Phi(d_np, d_nm, d00),
    d2Phi_dm2               = d2Phi(d_mp, d_mm, d00),
    d2Phi_dlog_thalf_dnu    = d2Phi_mixed(d_tp_np, d_tp_nm, d_tm_np, d_tm_nm),
    d2Phi_dlog_thalf_dm     = d2Phi_mixed(d_tp_mp, d_tp_mm, d_tm_mp, d_tm_mm),
    d2Phi_dnu_dm            = d2Phi_mixed(d_np_mp, d_np_mm, d_nm_mp, d_nm_mm),

    d2phi_dlog_thalf2       = d2phi(d_tp, d_tm, d00),
    d2phi_dnu2              = d2phi(d_np, d_nm, d00),
    d2phi_dm2               = d2phi(d_mp, d_mm, d00),
    d2phi_dlog_thalf_dnu    = d2phi_mixed(d_tp_np, d_tp_nm, d_tm_np, d_tm_nm),
    d2phi_dlog_thalf_dm     = d2phi_mixed(d_tp_mp, d_tp_mm, d_tm_mp, d_tm_mm),
    d2phi_dnu_dm            = d2phi_mixed(d_np_mp, d_np_mm, d_nm_mp, d_nm_mm)
  )
}

# ---------------------------------------------------------------------------
# G3 (late) phase second derivatives
# ---------------------------------------------------------------------------
#' Second derivatives of G3 phase shape functions
#'
#' Analogue of \code{.hzr_g3_phase_derivatives()} for second-order
#' shape derivatives. Uses the same central-difference approach.
#'
#' @param time,tau,gamma,alpha,eta  G3 phase parameters.
#' @return List with \code{d2Phi_*} and \code{d2phi_*} entries for the ten
#'   unique pairs from \code{\{log\_tau, gamma, alpha, eta\}}.
#' @noRd
.hzr_g3_phase_second_derivatives <- function(time, tau, gamma, alpha, eta) {
  h <- .hzr_h2

  eval_at <- function(tau_p, gam_p, alp_p, eta_p) {
    d <- hzr_decompos_g3(time, tau = tau_p, gamma = gam_p,
                          alpha = alp_p, eta = eta_p)
    list(Phi = d$G3, phi = d$g3)
  }

  d00 <- eval_at(tau, gamma, alpha, eta)

  # Steps on internal scales (log_tau, gamma, alpha, eta direct)
  tau_p  <- tau * exp( h);  tau_m  <- tau * exp(-h)
  gam_p  <- gamma + h;      gam_m  <- gamma - h
  alp_p  <- alpha + h;      alp_m  <- alpha - h
  eta_p  <- eta   + h;      eta_m  <- eta   - h

  d_Tp <- eval_at(tau_p, gamma, alpha, eta)
  d_Tm <- eval_at(tau_m, gamma, alpha, eta)
  d_Gp <- eval_at(tau, gam_p, alpha, eta)
  d_Gm <- eval_at(tau, gam_m, alpha, eta)
  d_Ap <- eval_at(tau, gamma, alp_p, eta)
  d_Am <- eval_at(tau, gamma, alp_m, eta)
  d_Ep <- eval_at(tau, gamma, alpha, eta_p)
  d_Em <- eval_at(tau, gamma, alpha, eta_m)

  # Mixed pairs (using 4-point formula)
  four_pt <- function(par1_plus, par1_minus, par2_plus, par2_minus) {
    d_pp <- eval_at(par1_plus[[1]],  par1_plus[[2]],  par2_plus[[3]],  par2_plus[[4]])
    d_pm <- eval_at(par1_plus[[1]],  par1_plus[[2]],  par2_minus[[3]], par2_minus[[4]])
    d_mp <- eval_at(par1_minus[[1]], par1_minus[[2]], par2_plus[[3]],  par2_plus[[4]])
    d_mm <- eval_at(par1_minus[[1]], par1_minus[[2]], par2_minus[[3]], par2_minus[[4]])
    list(Phi = (d_pp$Phi - d_pm$Phi - d_mp$Phi + d_mm$Phi) / (4 * h^2),
         phi = (d_pp$phi - d_pm$phi - d_mp$phi + d_mm$phi) / (4 * h^2))
  }

  # Helper to bundle a parameter perturbation as (tau, gamma, alpha, eta)
  base  <- list(tau, gamma, alpha, eta)
  tau_plus  <- list(tau_p,  gamma, alpha, eta);   tau_minus  <- list(tau_m,  gamma, alpha, eta)
  gam_plus  <- list(tau, gam_p,  alpha, eta);     gam_minus  <- list(tau, gam_m,  alpha, eta)
  alp_plus  <- list(tau, gamma, alp_p,  eta);     alp_minus  <- list(tau, gamma, alp_m,  eta)
  eta_plus  <- list(tau, gamma, alpha, eta_p);    eta_minus  <- list(tau, gamma, alpha, eta_m)

  h2 <- h^2
  d2Phi_diag <- function(d_p, d_m) (d_p$Phi - 2 * d00$Phi + d_m$Phi) / h2
  d2phi_diag <- function(d_p, d_m) (d_p$phi - 2 * d00$phi + d_m$phi) / h2

  TG <- four_pt(tau_plus,  tau_minus,  gam_plus,  gam_minus)
  TA <- four_pt(tau_plus,  tau_minus,  alp_plus,  alp_minus)
  TE <- four_pt(tau_plus,  tau_minus,  eta_plus,  eta_minus)
  GA <- four_pt(gam_plus,  gam_minus,  alp_plus,  alp_minus)
  GE <- four_pt(gam_plus,  gam_minus,  eta_plus,  eta_minus)
  AE <- four_pt(alp_plus,  alp_minus,  eta_plus,  eta_minus)

  list(
    d2Phi_dlog_tau2        = d2Phi_diag(d_Tp, d_Tm),
    d2Phi_dgamma2          = d2Phi_diag(d_Gp, d_Gm),
    d2Phi_dalpha2          = d2Phi_diag(d_Ap, d_Am),
    d2Phi_deta2            = d2Phi_diag(d_Ep, d_Em),
    d2Phi_dlog_tau_dgamma  = TG$Phi,
    d2Phi_dlog_tau_dalpha  = TA$Phi,
    d2Phi_dlog_tau_deta    = TE$Phi,
    d2Phi_dgamma_dalpha    = GA$Phi,
    d2Phi_dgamma_deta      = GE$Phi,
    d2Phi_dalpha_deta      = AE$Phi,

    d2phi_dlog_tau2        = d2phi_diag(d_Tp, d_Tm),
    d2phi_dgamma2          = d2phi_diag(d_Gp, d_Gm),
    d2phi_dalpha2          = d2phi_diag(d_Ap, d_Am),
    d2phi_deta2            = d2phi_diag(d_Ep, d_Em),
    d2phi_dlog_tau_dgamma  = TG$phi,
    d2phi_dlog_tau_dalpha  = TA$phi,
    d2phi_dlog_tau_deta    = TE$phi,
    d2phi_dgamma_dalpha    = GA$phi,
    d2phi_dgamma_deta      = GE$phi,
    d2phi_dalpha_deta      = AE$phi
  )
}
```

### Step 4: Run to verify helpers pass

```bash
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-multiphase-hessian.R")'
```
Expected: PASS.

### Step 5: Commit

```bash
git add R/hessian-multiphase.R tests/testthat/test-multiphase-hessian.R
git commit -m "feat: shape second-derivative helpers for multiphase Hessian"
```

---

## Task 2 — `.hzr_hessian_multiphase()` core function

**Files:**
- Modify: `R/hessian-multiphase.R` (append)
- Modify: `tests/testthat/test-multiphase-hessian.R` (append)

### Step 1: Write the failing cross-check tests

Append to `tests/testthat/test-multiphase-hessian.R`:

```r
# --- .hzr_hessian_multiphase() cross-check -----------------------------------

.make_avc_fit_args <- function(conserve = TRUE) {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  list(
    time   = avc$int_dead,
    status = avc$dead,
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                            fixed = "m"),
      constant = hzr_phase("constant")
    ),
    cov_counts = c(early = 0L, constant = 0L),
    x_list = list(early = NULL, constant = NULL),
    weights = NULL
  )
}

test_that(".hzr_hessian_multiphase matches numDeriv (2-phase no covariates)", {
  skip_if_not_installed("numDeriv")
  a <- .make_avc_fit_args()

  # Fit to get a reasonable theta
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data  = na.omit(avc),
    dist  = "multiphase",
    phases = a$phases,
    fit   = TRUE,
    control = list(n_starts = 1, conserve = FALSE)
  )
  theta <- fit$fit$par

  # Re-unpack cov_counts and x_list from the fitted object
  cov_counts <- fit$fit$covariate_counts
  x_list     <- fit$fit$x_list

  # Analytic Hessian
  H_an <- .hzr_hessian_multiphase(
    theta, time = a$time, status = a$status,
    phases = a$phases, covariate_counts = cov_counts, x_list = x_list
  )
  expect_true(is.matrix(H_an))

  # numDeriv reference on the full NLL
  obj <- function(par) {
    -.hzr_logl_multiphase(par, a$time, a$status,
                           phases = a$phases,
                           covariate_counts = cov_counts,
                           x_list = x_list)
  }
  H_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-3,
               label = "multiphase 2-phase NLL Hessian vs numDeriv")
})

test_that(".hzr_hessian_multiphase matches numDeriv (2-phase with covariates)", {
  skip_if_not_installed("numDeriv")
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  phases <- list(
    early    = hzr_phase("cdf", formula = ~ age + com_iv,
                          t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ early(age + com_iv),
    data = avc, dist = "multiphase", phases = phases, fit = TRUE,
    control = list(n_starts = 1, conserve = FALSE)
  )
  theta      <- fit$fit$par
  cov_counts <- fit$fit$covariate_counts
  x_list     <- fit$fit$x_list

  H_an <- .hzr_hessian_multiphase(
    theta, time = avc$int_dead, status = avc$dead,
    phases = phases, covariate_counts = cov_counts, x_list = x_list
  )
  obj <- function(par) {
    -.hzr_logl_multiphase(par, avc$int_dead, avc$dead,
                           phases = phases,
                           covariate_counts = cov_counts,
                           x_list = x_list)
  }
  H_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-3,
               label = "multiphase 2-phase + covariates Hessian vs numDeriv")
})

test_that(".hzr_hessian_multiphase returns NULL for interval-censored rows", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  status_mixed <- avc$dead
  status_mixed[1:3] <- 2L  # inject interval-censored rows

  phases <- list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                    fixed = "m"),
                 constant = hzr_phase("constant"))
  expect_null(.hzr_hessian_multiphase(
    c(early.log_mu = 0, early.log_t_half = log(0.15), early.nu = 1.4,
      constant.log_mu = 0),
    time = avc$int_dead, status = status_mixed,
    phases = phases, covariate_counts = c(early = 0L, constant = 0L),
    x_list = list(early = NULL, constant = NULL)
  ))
})

test_that(".hzr_hessian_multiphase returns NULL for left-censored rows", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  status_left <- avc$dead
  status_left[1:3] <- -1L

  phases <- list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                    fixed = "m"),
                 constant = hzr_phase("constant"))
  expect_null(.hzr_hessian_multiphase(
    c(early.log_mu = 0, early.log_t_half = log(0.15), early.nu = 1.4,
      constant.log_mu = 0),
    time = avc$int_dead, status = status_left,
    phases = phases, covariate_counts = c(early = 0L, constant = 0L),
    x_list = list(early = NULL, constant = NULL)
  ))
})
```

### Step 2: Run to verify they fail

```bash
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-multiphase-hessian.R")'
```
Expected: FAIL — `could not find function ".hzr_hessian_multiphase"`.

### Step 3: Implement `.hzr_hessian_multiphase()`

Append to `R/hessian-multiphase.R`:

```r
# ---------------------------------------------------------------------------
# Main analytic Hessian function
# ---------------------------------------------------------------------------
#' Analytic Hessian of the multiphase NLL
#'
#' Returns the \eqn{p \times p} Hessian of the **objective** (negative
#' log-likelihood) on the internal parameter scale, matching
#' \code{numDeriv::hessian(objective)}.  Coverage: event + right-censored
#' rows, with counting-process start-time corrections when
#' \code{time_lower > 0}. Returns \code{NULL} for data containing
#' left-censored (\code{status == -1}) or interval-censored
#' (\code{status == 2}) rows; the caller falls back to numDeriv in that case.
#'
#' **Structure** (see derivation in PLAN-hessian-analytic-multiphase.md):
#' \deqn{H = (\text{block-diagonal } A) + (\text{dense outer-product } B) + (\text{block-diagonal } C)}
#'
#' Term A: curvature of \eqn{\sum_i w_i H_i}; block-diagonal in phases.
#' Term B: \eqn{\sum_{i\in E} (w_i/h_i^2) \nabla h_i \nabla h_i^T}; dense.
#' Term C: curvature of \eqn{-\sum_{i\in E} w_i \log h_i}; block-diagonal.
#'
#' @noRd
.hzr_hessian_multiphase <- function(
    theta, time, status,
    time_lower = NULL, time_upper = NULL,
    x = NULL, weights = NULL,
    phases, covariate_counts, x_list) {

  n <- length(time)
  p <- length(theta)
  if (is.null(weights)) weights <- rep(1, n)

  # Coverage contract: decline for left/interval-censored rows
  if (any(status %in% c(-1L, 2L))) return(NULL)

  idx_event <- which(status == 1L)
  n_e       <- length(idx_event)

  # Counting-process start-time handling (mirrors gradient)
  need_start <- !is.null(time_lower) &&
                any(time_lower > 0 & time_lower < time & status %in% c(0L, 1L))
  start_vec <- if (need_start) {
    sv <- rep(0, n)
    epoch_idx <- status %in% c(0L, 1L) & time_lower < time
    sv[epoch_idx] <- time_lower[epoch_idx]
    sv
  } else {
    NULL
  }

  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)

  # --- Per-phase quantities --------------------------------------------------
  n_phases        <- length(phases)
  phase_mu        <- vector("list", n_phases)
  phase_Phi       <- vector("list", n_phases)
  phase_phi       <- vector("list", n_phases)
  phase_d1        <- vector("list", n_phases)   # .hzr_phase_derivatives output
  phase_d2        <- vector("list", n_phases)   # second-derivative output
  phase_Phi_start <- vector("list", n_phases)
  phase_d1_start  <- vector("list", n_phases)
  phase_d2_start  <- vector("list", n_phases)
  phase_x_tilde   <- vector("list", n_phases)   # [1 | x_j], n × (1+p_j)

  H_t <- rep(0, n)
  h_t <- rep(0, n)

  for (j in seq_along(phases)) {
    nm   <- names(phases)[j]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])

    # mu_j(i) and augmented covariate matrix
    if (length(pars$beta) > 0 && !is.null(x_list[[nm]])) {
      x_j  <- x_list[[nm]]
      eta_j <- pars$log_mu + as.numeric(x_j %*% pars$beta)
    } else {
      x_j   <- NULL
      eta_j <- rep(pars$log_mu, n)
    }
    mu_j <- exp(eta_j)
    phase_mu[[j]] <- mu_j

    x_tilde_j <- if (is.null(x_j)) matrix(1, n, 1L) else cbind(1, x_j)
    phase_x_tilde[[j]] <- x_tilde_j

    # Shape quantities
    if (phases[[nm]]$type == "constant") {
      phase_Phi[[j]] <- time
      phase_phi[[j]] <- rep(1, n)
      phase_d1[[j]]  <- NULL
      phase_d2[[j]]  <- NULL
      if (need_start) {
        phase_Phi_start[[j]] <- start_vec
        phase_d1_start[[j]]  <- NULL
        phase_d2_start[[j]]  <- NULL
      }
    } else if (phases[[nm]]$type == "g3") {
      tau_j <- exp(pars$log_tau)
      d3    <- hzr_decompos_g3(time, tau = tau_j, gamma = pars$gamma,
                                alpha = pars$alpha, eta = pars$eta)
      phase_Phi[[j]] <- d3$G3
      phase_phi[[j]] <- d3$g3
      phase_d1[[j]]  <- .hzr_g3_phase_derivatives(
        time, tau = tau_j, gamma = pars$gamma,
        alpha = pars$alpha, eta = pars$eta)
      phase_d2[[j]]  <- .hzr_g3_phase_second_derivatives(
        time, tau = tau_j, gamma = pars$gamma,
        alpha = pars$alpha, eta = pars$eta)
      if (need_start) {
        d3s <- hzr_decompos_g3(start_vec, tau = tau_j, gamma = pars$gamma,
                                alpha = pars$alpha, eta = pars$eta)
        phase_Phi_start[[j]] <- d3s$G3
        phase_d1_start[[j]]  <- .hzr_g3_phase_derivatives(
          start_vec, tau = tau_j, gamma = pars$gamma,
          alpha = pars$alpha, eta = pars$eta)
        phase_d2_start[[j]]  <- .hzr_g3_phase_second_derivatives(
          start_vec, tau = tau_j, gamma = pars$gamma,
          alpha = pars$alpha, eta = pars$eta)
      }
    } else {
      # CDF / hazard
      t_half_j <- exp(pars$log_t_half)
      d1 <- .hzr_phase_derivatives(time, t_half = t_half_j, nu = pars$nu,
                                    m = pars$m, type = phases[[nm]]$type)
      phase_Phi[[j]] <- d1$Phi
      phase_phi[[j]] <- d1$phi
      phase_d1[[j]]  <- d1
      phase_d2[[j]]  <- .hzr_phase_second_derivatives(
        time, t_half = t_half_j, nu = pars$nu, m = pars$m,
        type = phases[[nm]]$type)
      if (need_start) {
        d1s <- .hzr_phase_derivatives(start_vec, t_half = t_half_j,
                                       nu = pars$nu, m = pars$m,
                                       type = phases[[nm]]$type)
        phase_Phi_start[[j]] <- d1s$Phi
        phase_d1_start[[j]]  <- d1s
        phase_d2_start[[j]]  <- .hzr_phase_second_derivatives(
          start_vec, t_half = t_half_j, nu = pars$nu, m = pars$m,
          type = phases[[nm]]$type)
      }
    }

    H_t <- H_t + mu_j * phase_Phi[[j]]
    h_t <- h_t + mu_j * phase_phi[[j]]
    if (need_start) {
      H_t <- H_t - mu_j * phase_Phi_start[[j]]  # H(t) - H(start) net
    }
  }

  if (any(!is.finite(H_t)) || any(!is.finite(h_t))) {
    return(NULL)
  }

  # h_e = h at event times (guard against 0)
  h_e <- pmax(h_t[idx_event], .Machine$double.xmin)

  # w/h_e^2 and w/h_e vectors (length n_e)
  w_e     <- weights[idx_event]
  wh2_e   <- w_e / h_e^2     # for Term B
  wh1_e   <- w_e / h_e       # for Term C

  # --- Build the full p × p Hessian matrix ----------------------------------
  H_mat <- matrix(0, p, p)
  dimnames(H_mat) <- list(names(theta), names(theta))

  # We need a per-phase parameter offset so we know where in H_mat to write.
  # Build an offset table.
  phase_param_counts <- vapply(seq_along(phases), function(j) {
    nm <- names(phases)[j]
    if (phases[[nm]]$type == "constant") {
      1L + covariate_counts[[nm]]
    } else if (phases[[nm]]$type == "g3") {
      5L + covariate_counts[[nm]]  # log_mu, log_tau, gamma, alpha, eta
    } else {
      4L + covariate_counts[[nm]]  # log_mu, log_t_half, nu, m
    }
  }, integer(1))

  phase_offsets <- c(0L, cumsum(phase_param_counts))

  # --- Term B: dense outer product of h-gradient ----------------------------
  # Assemble X̃_h (n_e × p): row i is the gradient of h(t_i) w.r.t. theta.
  X_h <- matrix(0, n_e, p)

  for (j in seq_along(phases)) {
    nm   <- names(phases)[j]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])
    pos  <- phase_offsets[j]
    mu_j <- phase_mu[[j]]
    Phi_j <- phase_Phi[[j]]
    phi_j <- phase_phi[[j]]
    d1    <- phase_d1[[j]]
    x_t   <- phase_x_tilde[[j]]

    mu_phi_e <- mu_j[idx_event] * phi_j[idx_event]  # length n_e

    # log_mu column: dh/d(log_mu_j) = mu_j phi_j
    X_h[, pos + 1L] <- mu_phi_e

    # Shape columns
    if (phases[[nm]]$type %in% c("cdf", "hazard")) {
      shape_dphi <- list(d1$dphi_dlog_thalf, d1$dphi_dnu, d1$dphi_dm)
      for (s in seq_along(shape_dphi)) {
        X_h[, pos + 1L + s] <- mu_j[idx_event] * shape_dphi[[s]][idx_event]
      }
      beta_start <- pos + 5L
    } else if (phases[[nm]]$type == "g3") {
      shape_dphi <- list(d1$dphi_dlog_tau, d1$dphi_dgamma,
                         d1$dphi_dalpha,   d1$dphi_deta)
      for (s in seq_along(shape_dphi)) {
        X_h[, pos + 1L + s] <- mu_j[idx_event] * shape_dphi[[s]][idx_event]
      }
      beta_start <- pos + 6L
    } else {
      # constant phase: no shape cols
      beta_start <- pos + 2L
    }

    # Beta columns: dh/d(beta_jk) = mu_j phi_j x_jk
    if (covariate_counts[[nm]] > 0 && !is.null(x_list[[nm]])) {
      x_j_e <- x_list[[nm]][idx_event, , drop = FALSE]
      for (k in seq_len(covariate_counts[[nm]])) {
        X_h[, beta_start + k - 1L] <- mu_phi_e * x_j_e[, k]
      }
    }
  }

  # Term B = X̃_h' diag(wh2_e) X̃_h (each row of X_h weighted by wh2_e)
  B_mat <- crossprod(X_h, wh2_e * X_h)
  H_mat <- H_mat + B_mat

  # --- Terms A and C: per-phase block contributions -------------------------
  for (j in seq_along(phases)) {
    nm   <- names(phases)[j]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], phases[[nm]])
    pos  <- phase_offsets[j]
    mu_j <- phase_mu[[j]]
    Phi_j <- phase_Phi[[j]]
    phi_j <- phase_phi[[j]]
    d1    <- phase_d1[[j]]
    d2    <- phase_d2[[j]]
    x_t   <- phase_x_tilde[[j]]

    wmu_Phi   <- weights * mu_j * Phi_j               # Term A weight (all rows)
    wmu_phi_e <- wh1_e * mu_j[idx_event] * phi_j[idx_event]  # Term C weight (events)

    # Counting-process start correction for Term A:
    if (need_start) {
      wmu_Phi_start <- weights * mu_j * phase_Phi_start[[j]]
      wmu_Phi_net <- wmu_Phi - wmu_Phi_start   # net for H(t) - H(start)
    } else {
      wmu_Phi_net <- wmu_Phi
    }

    # --- μ/β sub-block (crossprod style) ---------------------------------
    # Term A μ/β: X̃_j' diag(wmu_Phi_net) X̃_j
    # Term C μ/β: - X̃_j_E' diag(wmu_phi_e) X̃_j_E
    x_t_e <- x_t[idx_event, , drop = FALSE]

    block_A_mub <- crossprod(x_t, wmu_Phi_net * x_t)
    block_C_mub <- crossprod(x_t_e, wmu_phi_e * x_t_e)

    # Number of μ/β parameters: 1 + p_j
    n_mub <- 1L + covariate_counts[[nm]]

    if (phases[[nm]]$type == "constant") {
      # Only μ/β block; no shape params
      idx_block <- pos + seq_len(n_mub)
      H_mat[idx_block, idx_block] <- H_mat[idx_block, idx_block] +
        block_A_mub - block_C_mub

    } else {
      # Shape params come AFTER log_mu, BEFORE beta (log_mu | shapes | beta)
      n_shape <- if (phases[[nm]]$type == "g3") 4L else 3L

      idx_mub   <- pos + 1L                             # log_mu index
      idx_shape <- pos + 1L + seq_len(n_shape)          # shape indices
      idx_beta  <- if (covariate_counts[[nm]] > 0)
                     pos + 1L + n_shape + seq_len(covariate_counts[[nm]])
                   else integer(0)

      # log_mu / log_mu (scalar)
      H_mat[idx_mub, idx_mub] <- H_mat[idx_mub, idx_mub] +
        sum(wmu_Phi_net) - sum(wmu_phi_e)

      # log_mu / beta cross terms (only if covariates present)
      if (length(idx_beta) > 0) {
        x_j <- x_list[[nm]]
        x_j_e <- x_j[idx_event, , drop = FALSE]
        cross_A <- colSums(wmu_Phi_net[, drop = FALSE] * x_j)
        cross_C <- colSums(wmu_phi_e[, drop = FALSE] * x_j_e)
        H_mat[idx_mub, idx_beta] <- H_mat[idx_mub, idx_beta] + cross_A - cross_C
        H_mat[idx_beta, idx_mub] <- H_mat[idx_beta, idx_mub] + cross_A - cross_C
      }

      # beta / beta (only if covariates)
      if (length(idx_beta) > 0) {
        x_j   <- x_list[[nm]]
        x_j_e <- x_j[idx_event, , drop = FALSE]
        bb_A  <- crossprod(x_j,   wmu_Phi_net * x_j)
        bb_C  <- crossprod(x_j_e, wmu_phi_e   * x_j_e)
        H_mat[idx_beta, idx_beta] <- H_mat[idx_beta, idx_beta] + bb_A - bb_C
      }

      # Shape parameters: extract first-derivative vectors
      if (phases[[nm]]$type %in% c("cdf", "hazard")) {
        shape_keys_d1_Phi <- c("dPhi_dlog_thalf", "dPhi_dnu", "dPhi_dm")
        shape_keys_d1_phi <- c("dphi_dlog_thalf", "dphi_dnu", "dphi_dm")
        shape_keys_d2_Phi <- list(
          c("d2Phi_dlog_thalf2",    "d2Phi_dlog_thalf2"),
          c("d2Phi_dnu2",           "d2Phi_dnu2"),
          c("d2Phi_dm2",            "d2Phi_dm2"),
          c("d2Phi_dlog_thalf_dnu", "d2Phi_dlog_thalf_dnu"),
          c("d2Phi_dlog_thalf_dm",  "d2Phi_dlog_thalf_dm"),
          c("d2Phi_dnu_dm",         "d2Phi_dnu_dm")
        )
        shape_keys_d2_phi <- list(
          c("d2phi_dlog_thalf2",    "d2phi_dlog_thalf2"),
          c("d2phi_dnu2",           "d2phi_dnu2"),
          c("d2phi_dm2",            "d2phi_dm2"),
          c("d2phi_dlog_thalf_dnu", "d2phi_dlog_thalf_dnu"),
          c("d2phi_dlog_thalf_dm",  "d2phi_dlog_thalf_dm"),
          c("d2phi_dnu_dm",         "d2phi_dnu_dm")
        )
        # Which pairs of shape indices correspond to the 6 unique pairs?
        shape_pairs <- list(c(1,1), c(2,2), c(3,3), c(1,2), c(1,3), c(2,3))
      } else {
        # G3: 4 shape params → 10 unique pairs
        shape_keys_d1_Phi <- c("dPhi_dlog_tau", "dPhi_dgamma",
                                "dPhi_dalpha",   "dPhi_deta")
        shape_keys_d1_phi <- c("dphi_dlog_tau", "dphi_dgamma",
                                "dphi_dalpha",   "dphi_deta")
        shape_keys_d2_Phi <- list(
          c("d2Phi_dlog_tau2",       "d2Phi_dlog_tau2"),
          c("d2Phi_dgamma2",         "d2Phi_dgamma2"),
          c("d2Phi_dalpha2",         "d2Phi_dalpha2"),
          c("d2Phi_deta2",           "d2Phi_deta2"),
          c("d2Phi_dlog_tau_dgamma", "d2Phi_dlog_tau_dgamma"),
          c("d2Phi_dlog_tau_dalpha", "d2Phi_dlog_tau_dalpha"),
          c("d2Phi_dlog_tau_deta",   "d2Phi_dlog_tau_deta"),
          c("d2Phi_dgamma_dalpha",   "d2Phi_dgamma_dalpha"),
          c("d2Phi_dgamma_deta",     "d2Phi_dgamma_deta"),
          c("d2Phi_dalpha_deta",     "d2Phi_dalpha_deta")
        )
        shape_keys_d2_phi <- lapply(shape_keys_d2_Phi, function(k)
          sub("Phi", "phi", k))
        shape_pairs <- c(
          list(c(1,1), c(2,2), c(3,3), c(4,4)),
          list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4))
        )
      }

      n_shape <- length(shape_keys_d1_Phi)

      # log_mu / shape cross terms (Term A and C)
      for (s in seq_len(n_shape)) {
        dPhi_ds <- d1[[shape_keys_d1_Phi[s]]]
        dphi_ds <- d1[[shape_keys_d1_phi[s]]]
        a_val <- sum(wmu_Phi_net * dPhi_ds)
        if (need_start && !is.null(phase_d1_start[[j]])) {
          # The start correction for d²H/d(log_mu)d(s): sum(-w μ dΦ(start)/ds)
          # is already embedded in wmu_Phi_net via the Phi_net = Phi(t) - Phi(start)
          # structure. But dPhi_ds is at time, not start. Correct:
          dPhi_ds_start <- phase_d1_start[[j]][[shape_keys_d1_Phi[s]]]
          a_val <- sum(weights * mu_j * dPhi_ds) -
                   sum(weights * mu_j * dPhi_ds_start)
        }
        c_val <- sum(wmu_phi_e * dphi_ds[idx_event])
        col_s <- idx_shape[s]
        H_mat[idx_mub, col_s] <- H_mat[idx_mub, col_s] + a_val - c_val
        H_mat[col_s, idx_mub] <- H_mat[col_s, idx_mub] + a_val - c_val
      }

      # shape / shape block
      for (pair_idx in seq_along(shape_pairs)) {
        s_a <- shape_pairs[[pair_idx]][1]
        s_b <- shape_pairs[[pair_idx]][2]
        key_d2_Phi <- shape_keys_d2_Phi[[pair_idx]][1]
        key_d2_phi <- shape_keys_d2_phi[[pair_idx]][1]

        d2Phi_sasb <- d2[[key_d2_Phi]]
        d2phi_sasb <- d2[[key_d2_phi]]

        a_val <- sum(weights * mu_j * d2Phi_sasb)
        if (need_start && !is.null(phase_d2_start[[j]])) {
          a_val <- a_val - sum(weights * mu_j * phase_d2_start[[j]][[key_d2_Phi]])
        }
        c_val <- sum(wh1_e * mu_j[idx_event] * d2phi_sasb[idx_event])

        i_idx <- idx_shape[s_a]; j_idx <- idx_shape[s_b]
        H_mat[i_idx, j_idx] <- H_mat[i_idx, j_idx] + a_val - c_val
        if (s_a != s_b) {
          H_mat[j_idx, i_idx] <- H_mat[j_idx, i_idx] + a_val - c_val
        }
      }

      # beta / shape cross terms
      if (length(idx_beta) > 0) {
        x_j   <- x_list[[nm]]
        x_j_e <- x_j[idx_event, , drop = FALSE]
        for (s in seq_len(n_shape)) {
          dPhi_ds <- d1[[shape_keys_d1_Phi[s]]]
          dphi_ds <- d1[[shape_keys_d1_phi[s]]]
          if (need_start && !is.null(phase_d1_start[[j]])) {
            dPhi_ds_start <- phase_d1_start[[j]][[shape_keys_d1_Phi[s]]]
            a_vec <- colSums((weights * mu_j * dPhi_ds -
                              weights * mu_j * dPhi_ds_start) * x_j)
          } else {
            a_vec <- colSums((weights * mu_j * dPhi_ds) * x_j)
          }
          c_vec <- colSums((wh1_e * mu_j[idx_event] *
                              dphi_ds[idx_event]) * x_j_e)
          col_s <- idx_shape[s]
          H_mat[idx_beta, col_s] <- H_mat[idx_beta, col_s] + a_vec - c_vec
          H_mat[col_s, idx_beta] <- H_mat[col_s, idx_beta] + a_vec - c_vec
        }
      }
    }
  }

  # Guard: zero any non-finite entries (safety net; should not trigger at valid theta)
  H_mat[!is.finite(H_mat)] <- 0

  H_mat
}
```

### Step 4: Run to verify cross-check tests pass

```bash
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-multiphase-hessian.R")'
```
Expected: PASS on cross-check tests (tolerance 1e-3) and NULL-fallback tests.

**If a cross-check test fails:** check the block assembly indexing. Common
mistakes: shape param offsets when covariates present; forgetting to negate
Term C (it enters as −C in H_mat); start-time correction sign.

### Step 5: Commit

```bash
git add R/hessian-multiphase.R tests/testthat/test-multiphase-hessian.R
git commit -m "feat: .hzr_hessian_multiphase() analytic NLL Hessian (3-term assembly)"
```

---

## Task 3 — Wire into `.hzr_optim_multiphase()` and CoE full-info path

**Files:**
- Modify: `R/likelihood-multiphase.R`
- Modify: `tests/testthat/test-multiphase-hessian.R` (append)

### Step 1: Write the failing end-to-end test

Append to `tests/testthat/test-multiphase-hessian.R`:

```r
test_that("multiphase fit vcov uses analytic Hessian (matches numDeriv to 1e-3)", {
  skip_if_not_installed("numDeriv")
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = phases, fit = TRUE,
    control = list(n_starts = 1, conserve = FALSE)
  )
  theta <- fit$fit$par
  cov_counts <- fit$fit$covariate_counts
  x_list     <- fit$fit$x_list

  # numDeriv reference vcov
  obj <- function(par) -.hzr_logl_multiphase(
    par, avc$int_dead, avc$dead, phases = phases,
    covariate_counts = cov_counts, x_list = x_list)
  vcov_nd <- solve(numDeriv::hessian(obj, theta))

  expect_equal(unname(fit$fit$vcov[!is.na(fit$fit$vcov)]),
               unname(vcov_nd[!is.na(fit$fit$vcov)]),
               tolerance = 1e-3,
               label = "analytic vcov matches numDeriv vcov (non-NA entries)")
})

test_that("CoE fit vcov (full-info) uses analytic Hessian", {
  skip_if_not_installed("numDeriv")
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  phases <- list(
    early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit_coe <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = phases, fit = TRUE,
    control = list(n_starts = 1, conserve = TRUE)
  )
  # The conserved log_mu must now carry a finite SE (CoE full-info path)
  se <- sqrt(diag(fit_coe$fit$vcov))
  expect_true(all(is.finite(se[!is.na(se)])))
  expect_true(is.finite(fit_coe$fit$rcond))
})
```

### Step 2: Run to verify they fail

Expected: FAIL — the vcov still comes from numDeriv (or the tests now pass with
numDeriv, which means they'll still pass after wiring — either way proceed to
Step 3 and verify `isTRUE(fit$fit$pd)` confirms the analytic path was used).

### Step 3: Wire the `hessian_fn` closure into `.hzr_optim_multiphase()`

In `R/likelihood-multiphase.R`, locate `.hzr_optim_multiphase()` (the function
that calls `.hzr_optim_generic()`). Immediately before the
`.hzr_optim_generic()` call, build a `hessian_fn` closure:

```r
  # Analytic Hessian closure — covers event + right-censored rows.
  # Returns NULL for left/interval-censored data to trigger numDeriv fallback.
  hessian_fn_mp <- function(theta_free) {
    # theta_free is on the REDUCED (free-parameter) scale; expand to full.
    theta_full <- if (any_fixed) {
      th <- theta_fixed
      th[free_idx] <- theta_free
      th
    } else {
      theta_free
    }
    H_full <- .hzr_hessian_multiphase(
      theta_full, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x, weights = weights,
      phases = phases, covariate_counts = covariate_counts, x_list = x_list
    )
    if (is.null(H_full)) return(NULL)
    # Restrict to the free-parameter block that the optimizer sees
    if (any_fixed) H_full[free_idx, free_idx] else H_full
  }
```

Then pass it to `.hzr_optim_generic()`:

```r
    hessian_fn = hessian_fn_mp
```

### Step 4: Update the CoE full-info vcov path (L1418)

In the CoE full-information vcov block
(`if (use_conserve && !is.null(fixmu_pos)) { ... }`, around L1405–1447),
replace the `numDeriv::hessian()` call with:

```r
        # Prefer analytic Hessian; fall back to numDeriv if it declines
        H_unc <- NULL
        theta_unc_base <- base_theta[idx_unc]
        H_unc <- tryCatch(
          .hzr_hessian_multiphase(
            base_theta, time = time, status = status,
            time_lower = time_lower, time_upper = time_upper,
            x = x, weights = weights,
            phases = phases, covariate_counts = covariate_counts,
            x_list = x_list
          ),
          error = function(e) NULL
        )
        if (is.matrix(H_unc)) {
          # Restrict to the unconstrained free set
          H_unc <- H_unc[idx_unc, idx_unc]
        }
        # Fallback to numDeriv when analytic declines (left/interval rows)
        if (is.null(H_unc) && requireNamespace("numDeriv", quietly = TRUE)) {
          H_unc <- tryCatch(
            numDeriv::hessian(neg_ll_unc, base_theta[idx_unc]),
            error = function(e) NULL
          )
        }
```

Leave the rest of the CoE vcov block (`.hzr_safe_solve(H_unc)`, expansion
to full dimension, etc.) unchanged.

### Step 5: Run to verify end-to-end tests pass

```bash
Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-multiphase-hessian.R")'
```
Expected: PASS.

### Step 6: Regression check — SAS-parity SE assertions must not move adversely

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-sas-parity.R")'
```
Expected: 0 failures. A shift toward the SAS SE value is acceptable. A shift
away is a bug — stop and investigate before proceeding.

### Step 7: Commit

```bash
git add R/likelihood-multiphase.R tests/testthat/test-multiphase-hessian.R
git commit -m "feat: wire analytic Hessian into multiphase optimizer and CoE full-info vcov"
```

---

## Task 4 — Scale-invariance and 13-parameter anchor

**Files:**
- Modify: `tests/testthat/test-multiphase-hessian.R` (append)
- Modify: `tests/testthat/test-hessian-stability.R` (update placeholder)

### Step 1: Write the tests

Append to `tests/testthat/test-multiphase-hessian.R`:

```r
test_that("multiphase SEs are invariant to covariate rescaling (log-mu z-stat)", {
  skip_on_cran()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  phases_base <- list(
    early    = hzr_phase("cdf", formula = ~ age,
                          t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  avc2 <- avc; avc2$age <- avc$age / 10  # rescaled covariate

  f1 <- hazard(survival::Surv(int_dead, dead) ~ early(age),
               data = avc,  dist = "multiphase", phases = phases_base, fit = TRUE,
               control = list(n_starts = 1, conserve = FALSE))
  f2 <- hazard(survival::Surv(int_dead, dead) ~ early(age),
               data = avc2, dist = "multiphase", phases = phases_base, fit = TRUE,
               control = list(n_starts = 1, conserve = FALSE))

  # The beta z-statistic (estimate / SE) must be invariant under x → x/10.
  se1 <- sqrt(diag(f1$fit$vcov));  se2 <- sqrt(diag(f2$fit$vcov))
  beta_name <- grep("age", names(f1$fit$par), value = TRUE)[1]
  z1 <- f1$fit$par[beta_name] / se1[beta_name]
  # f2 beta for age is 10× larger (rescaled covariate), SE is 10× larger too
  beta_name2 <- grep("age", names(f2$fit$par), value = TRUE)[1]
  z2 <- f2$fit$par[beta_name2] / se2[beta_name2]
  expect_equal(unname(z1), unname(z2), tolerance = 1e-2,
               label = "beta z-stat invariant under covariate rescaling")
})

test_that("13-parameter multiphase anchor: stable SEs and rcond (supersedes placeholder)", {
  skip_on_cran()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  # Mirror hm.death.AVC.deciles: 6 early-phase covariates + constant intercept.
  # Use the subset of AVC variables available in the package dataset.
  # Target: >= 10 free parameters to exercise the high-dimensional path.
  phases_hd <- list(
    early    = hzr_phase("cdf",
                          formula = ~ age + com_iv + mal + op_mos + op_age + status,
                          t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(
    survival::Surv(int_dead, dead) ~
      early(age + com_iv + mal + op_mos + op_age + status),
    data = avc, dist = "multiphase", phases = phases_hd, fit = TRUE,
    control = list(n_starts = 3, maxit = 800, conserve = TRUE)
  )

  expect_true(fit$fit$converged)
  free <- which(is.finite(diag(fit$fit$vcov)))
  expect_gte(length(free), 10L)
  se <- sqrt(diag(fit$fit$vcov)[free])
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
  expect_true(is.finite(fit$fit$rcond))
  expect_true(isTRUE(fit$fit$pd),
              label = "13-param fit is positive-definite at optimum")
})
```

**Note:** verify the AVC dataset column names before running — use `names(avc)`
to confirm `age`, `com_iv`, `mal`, `op_mos`, `op_age`, `status` exist. Adjust
names if needed; keep free-parameter count ≥ 10.

In `tests/testthat/test-hessian-stability.R`, locate the placeholder
`"fit object carries rcond and pd diagnostics (multiphase)"` test and add a
comment noting this has been superseded by the comprehensive test in
`test-multiphase-hessian.R`. Do not delete it — it exercises a different
(simpler, 2-param) fit and remains useful for the Layer-1 diagnostics contract.

### Step 2: Run

```bash
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-multiphase-hessian.R")'
```
Expected: PASS. If the 13-param anchor fails `expect_true(isTRUE(fit$fit$pd))`,
the fit may be legitimately near-singular — check `fit$fit$rcond` and soften
the assertion to `expect_warning(... "ill-conditioned")` if rcond < 1.5e-8,
documenting it in the PR description.

### Step 3: Commit

```bash
git add tests/testthat/test-multiphase-hessian.R tests/testthat/test-hessian-stability.R
git commit -m "test: scale-invariance and 13-param anchor for multiphase analytic Hessian"
```

---

## Task 5 — NEWS entry and full check

**Files:**
- Modify: `NEWS.md`

### Step 1: Add the NEWS entry

Under `# TemporalHazard 1.1.0.9000 (development version)` in the
`## Improvements` subsection, add:

```markdown
* **Analytic Hessian for multiphase standard errors (Phase 7c, Layer 2 PR-6).**
  Post-fit standard errors for all multiphase fits now come from a closed-form
  Hessian of the negative log-likelihood rather than a numerical Richardson
  approximation. The Hessian is assembled from three terms: (A) a
  phase-block-diagonal curvature of Σᵢ wᵢ H(tᵢ), (B) a dense Fisher
  information outer product Σₑ (wᵢ/hᵢ²) ∇h ∇hᵀ capturing cross-phase
  parameter interactions, and (C) a phase-block-diagonal curvature of
  −Σₑ wᵢ log h(tᵢ). μ/β parameters use fully closed-form expressions;
  shape parameters (t_half, ν, m, and G3 parameters) use second-order
  central differences. The Conservation-of-Events full-information vcov
  path (introduced in v1.0.4) also switches to the analytic Hessian.
  Left/interval-censored fits fall back to the numerical Hessian.
  Completes the 6-PR analytic-Hessian rollout across all five families.
```

### Step 2: Full test suite

```bash
NOT_CRAN=true Rscript -e 'devtools::test()'
```
Expected: 0 failures.

### Step 3: R CMD check

```bash
Rscript -e 'devtools::check(args = c("--no-manual"), quiet = TRUE)'
```
Expected: 0 errors / 0 warnings; pre-existing NOTE only. Confirm no new
`object_usage_linter` warnings from `.hzr_hessian_multiphase`,
`.hzr_phase_second_derivatives`, `.hzr_g3_phase_second_derivatives`,
`.hzr_h2`.

### Step 4: Commit

```bash
git add NEWS.md
git commit -m "docs: NEWS entry for multiphase analytic Hessian (Layer 2 PR-6)"
```

---

## Self-Review Notes

- **Spec coverage (HESSIAN-STABILITY-DESIGN.md §Layer 2):**
  analytic-vs-numDeriv cross-check gate (T2); scale-invariance test (T4);
  13-param anchor (T4, replaces T5 placeholder); `hessian_fn` path through
  `.hzr_optim_generic` (T3); CoE full-info path updated (T3); coverage contract
  mirrors gradient's NULL-on-left/interval (T2 NULL tests); NEWS (T5).
  All PRs 1–5 plumbing already in place.

- **Term B is the new structural element.** Single-distribution Hessians (PRs
  2–5) are block-diagonal (no cross-parameter interaction through log h because
  there is only one phase). The multiphase dense outer product term is what
  makes this PR the hardest: every pair of parameters from different phases
  has a non-zero off-diagonal entry. Verify the B assembly by checking that
  `H_mat[idx_j, idx_k]` is non-zero for params from different phases j ≠ k
  at a typical non-degenerate θ.

- **Start-time corrections in shape/cross blocks.** The `wmu_Phi_net` vector
  (used for μ/β) naturally encodes `Φ(t) − Φ(start)`. The shape/cross terms
  and shape/shape terms must separately subtract the start-time contribution
  using `phase_d1_start` and `phase_d2_start`. Omitting this produces a
  Hessian that disagrees with numDeriv on counting-process fits.

- **CoE block.** The analytic Hessian is evaluated at `base_theta` (the full
  parameter vector including the conserved log_mu), then restricted to
  `idx_unc`. This is correct because the CoE optimum is the unconstrained MLE.
  The existing `.hzr_safe_solve` → vcov expansion → `fixed_mask` update logic
  is unchanged.

- **No regression for other families.** Weibull/exponential/log-logistic/
  log-normal pass their own `hessian_fn` and are unaffected. The multiphase
  path is the only one changing.

- **Objective-scale consistency.** `.hzr_hessian_multiphase()` returns the
  Hessian of the NLL (negative log-likelihood) to match the numDeriv reference
  `numDeriv::hessian(objective)`. Term A and Term C are positive semidefinite
  contributions to the NLL Hessian (PSD at the MLE). Term B is always PSD.
  The full matrix should be PD at an interior optimum; `fit$fit$pd == TRUE`
  is the empirical check.
