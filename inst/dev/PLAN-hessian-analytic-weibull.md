# Hessian Stability — Layer 2 PR-3 (Weibull Analytic Hessian) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an analytic Hessian as the production standard-error input for the
Weibull distribution, reusing the `hessian_fn` hook added in PR-2.

**Architecture:** `.hzr_optim_weibull()` fits on an internal, unconstrained
reparameterization `phi = (alpha, psi, beta)` where `g = exp(psi)`,
`H(t) = exp(alpha + eta) * t^g`, and then applies a delta-method `J` to map the
variance-covariance matrix to the natural `(mu, nu, beta)` scale. The analytic
Hessian must therefore be the **objective (negative log-likelihood) Hessian on
the internal `phi` scale** — matching `numDeriv::hessian(objective)`, which is
what `.hzr_optim_generic` currently inverts. The existing delta-method transform
downstream is unchanged. Coverage mirrors the analytic gradient: event +
right-censored only; return `NULL` on any left/interval-censored row so the
optimizer falls back to numDeriv.

**Tech Stack:** R, testthat (edition 3), roxygen2, devtools. Base-R numerics;
cross-check tests use `numDeriv` (Suggests; guard with `skip_if_not_installed`).

**Scope:** PR-3 of the Layer-2 sequence in `inst/dev/HESSIAN-STABILITY-DESIGN.md`.
Weibull only. The `hessian_fn` plumbing already exists (PR-2). Prerequisite:
PR #58/#59 (Weibull event-hazard + weighted-gradient fixes) are merged to `dev`,
so `logl_internal` and the natural-scale `.hzr_logl_weibull` agree for
event+right data. Branch: `dev` → v1.1.0.

---

## Mathematical derivation (internal scale)

Internal parameters `theta = (alpha, psi, beta_1..beta_p)` (order: 1 = alpha,
2 = psi, 3.. = beta). `g = exp(psi)`. Per row (event delta=1, right-censored
delta=0):

    eta_i  = x_i . beta           (0 if no covariates)
    H_i    = exp(alpha + eta_i) * t_i^g
    Hs_i   = exp(alpha + eta_i) * start_i^g      (epoch entry; start_i = 0 except
             status in {0,1} with time_lower < time, then start_i = time_lower)
    A_i    = H_i - Hs_i
    L_i    = log(t_i) ;  Ls_i = log(start_i)  (Ls_i = 0 where start_i = 0)

Log-likelihood (event + right):

    LL = sum_i w_i * [ delta_i * (alpha + eta_i + psi + (g-1) L_i) - A_i ]

The **objective** is `-LL`; `.hzr_safe_solve` inverts the objective Hessian.

Derivatives of H_i: `dH/dalpha = H`, `dH/dbeta_j = H x_ij`, `dH/dpsi = g L_i H_i`
(and analogously for Hs with Ls). Writing the design with an intercept column
`Xtilde = [1 | X]` (n x q, q = 1 + p_cov; column 1 is the alpha intercept):

**Objective Hessian blocks** (with `wA = w*A`, `wD = w*(H*L - Hs*Ls)`):

- `[alpha,beta] x [alpha,beta]` block (q x q):  `+ Xtilde' diag(wA) Xtilde`
- `psi` cross `[alpha,beta]` (length q):         `+ g * Xtilde' wD`
- `psi,psi` scalar:
  `- g*( sum(w*delta*L) - sum(w*H*L) + sum(w*Hs*Ls) )
   + g^2*( sum(w*H*L^2) - sum(w*Hs*Ls^2) )`

Assemble into the p x p matrix (p = 2 + p_cov) at index order (alpha=1, psi=2,
beta=3..p): the Xtilde block occupies rows/cols `idx_q = c(1, 3:p)`, psi is at 2.

Sanity: the `[alpha,beta]` block is exactly the exponential pattern
(`Xtilde' diag(w*cumhaz) Xtilde`) with the net cumulative hazard `A` and an
extra `psi` (shape) row/column. The numDeriv cross-check (Task 1) is the gate
that validates the algebra at arbitrary `phi`.

---

## File Structure

- **Modify** `R/likelihood-weibull.R` — add `.hzr_hessian_weibull_internal()`;
  wire it into `.hzr_optim_weibull()` via a `hessian_fn` closure.
- **Test** `tests/testthat/test-weibull-hessian.R` (new) — analytic-vs-numDeriv
  cross-check on the internal objective (no-cov; cov+weights; counting-process
  start), NULL fallback on left/interval, end-to-end natural-scale vcov match,
  scale-invariance.
- **Modify** `NEWS.md`.

---

## Task 1: `.hzr_hessian_weibull_internal()` analytic Hessian

**Files:**
- Modify: `R/likelihood-weibull.R` (add function before `.hzr_optim_weibull`)
- Test: `tests/testthat/test-weibull-hessian.R` (create)

- [ ] **Step 1: Write the failing cross-check tests**

Create `tests/testthat/test-weibull-hessian.R`:

```r
# Reconstruct the optimizer's internal objective (negative log-likelihood on the
# (alpha, psi, beta) scale). Post-PR#58 the natural-scale .hzr_logl_weibull
# equals logl_internal for event+right data, so we can map phi -> (mu, nu, beta)
# and reuse it.
.wb_obj_internal <- function(phi, time, status, time_lower = NULL,
                             time_upper = NULL, x = NULL, weights = NULL) {
  alpha <- phi[1]
  g <- exp(phi[2])
  beta <- if (length(phi) > 2) phi[3:length(phi)] else numeric(0)
  mu <- exp(alpha / g)
  theta_nat <- c(mu, g, beta)
  -.hzr_logl_weibull(theta_nat, time, status, time_lower, time_upper,
                     x = x, weights = weights)
}

test_that(".hzr_hessian_weibull_internal matches numDeriv (no covariates)", {
  skip_if_not_installed("numDeriv")
  set.seed(21)
  n <- 300
  time <- rexp(n, 0.6) + 0.01
  status <- rbinom(n, 1, 0.75)
  phi <- c(alpha = 0.4, psi = log(1.5))
  obj <- function(p) .wb_obj_internal(p, time, status)
  H_an <- .hzr_hessian_weibull_internal(phi, time, status)
  H_nd <- numDeriv::hessian(obj, phi)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_weibull_internal matches numDeriv (covariates + weights)", {
  skip_if_not_installed("numDeriv")
  set.seed(22)
  n <- 400
  x <- cbind(a = rnorm(n), b = rbinom(n, 1, 0.4))
  eta <- as.numeric(x %*% c(0.3, -0.4))
  time <- rweibull(n, shape = 1.4, scale = exp(-eta / 1.4)) + 0.01
  status <- rbinom(n, 1, 0.8)
  w <- runif(n, 0.5, 2)
  phi <- c(alpha = 0.2, psi = log(1.4), 0.3, -0.4)
  obj <- function(p) .wb_obj_internal(p, time, status, x = x, weights = w)
  H_an <- .hzr_hessian_weibull_internal(phi, time, status, x = x, weights = w)
  H_nd <- numDeriv::hessian(obj, phi)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_weibull_internal matches numDeriv (counting-process start)", {
  skip_if_not_installed("numDeriv")
  set.seed(23)
  n <- 250
  start <- runif(n, 0, 0.5)
  time <- start + rexp(n, 0.7) + 0.05
  status <- rbinom(n, 1, 0.7)
  phi <- c(alpha = 0.3, psi = log(1.2))
  obj <- function(p) .wb_obj_internal(p, time, status, time_lower = start)
  H_an <- .hzr_hessian_weibull_internal(phi, time, status, time_lower = start)
  H_nd <- numDeriv::hessian(obj, phi)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_weibull_internal returns NULL for left/interval censoring", {
  set.seed(24)
  n <- 40
  time <- rexp(n, 0.5) + 0.01
  status <- rep(1L, n)
  status[1:5] <- 2L
  expect_null(.hzr_hessian_weibull_internal(
    c(alpha = 0, psi = 0), time, status,
    time_lower = time, time_upper = time + 0.5))
  status2 <- rep(1L, n)
  status2[1:5] <- -1L
  expect_null(.hzr_hessian_weibull_internal(
    c(alpha = 0, psi = 0), time, status2, time_upper = time + 0.5))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-weibull-hessian.R")'`
Expected: FAIL — `could not find function ".hzr_hessian_weibull_internal"`.

- [ ] **Step 3: Implement `.hzr_hessian_weibull_internal()`**

In `R/likelihood-weibull.R`, immediately before `.hzr_optim_weibull()`, add:

```r
#' Analytic Hessian of the Weibull negative log-likelihood (internal scale)
#'
#' Returns the Hessian of the objective (negative log-likelihood) on the internal
#' \code{phi = (alpha, psi, beta)} parameterization that \code{.hzr_optim_weibull}
#' optimizes, where \code{g = exp(psi)} and \code{H(t) = exp(alpha + eta) t^g}.
#' This matches \code{numDeriv::hessian(objective)} so it can be used directly by
#' \code{.hzr_optim_generic()}; the downstream delta-method transform to the
#' natural \code{(mu, nu, beta)} scale is applied by the caller.  Coverage mirrors
#' the analytic gradient: event + right-censored rows only; returns \code{NULL}
#' for any left/interval-censored row (\code{status \%in\% c(-1, 2)}) to request
#' the numerical-Hessian fallback.
#'
#' @noRd
.hzr_hessian_weibull_internal <- function(
    theta, time, status,
    time_lower = NULL, time_upper = NULL, x = NULL, weights = NULL) {

  n <- length(time)
  if (is.null(weights)) weights <- rep(1, n)

  # Coverage contract: closed form is event + right censored only.
  if (any(status %in% c(-1, 2))) {
    return(NULL)
  }

  alpha <- theta[1]
  psi <- theta[2]
  g <- exp(psi)
  p <- length(theta)
  p_cov <- p - 2L
  beta <- if (p_cov > 0) theta[3:p] else numeric(0)
  eta <- if (p_cov > 0 && !is.null(x)) as.numeric(x %*% beta) else rep(0, n)

  # Cumulative hazards and counting-process entry-time (epochs only).
  base <- exp(alpha + eta)
  cumhaz <- base * (time ^ g)
  start_vec <- rep(0, n)
  if (!is.null(time_lower)) {
    epoch_idx <- status %in% c(0L, 1L) & time_lower < time
    start_vec[epoch_idx] <- time_lower[epoch_idx]
  }
  cumhaz_start <- base * (start_vec ^ g)

  log_t <- log(time)
  log_ts <- ifelse(start_vec > 0, log(start_vec), 0)

  w <- weights
  delta <- as.numeric(status == 1)
  net <- cumhaz - cumhaz_start                       # A_i

  # Weighted aggregate building blocks.
  wA   <- w * net
  wD   <- w * (cumhaz * log_t - cumhaz_start * log_ts)
  # psi,psi scalar.
  hpp <- -g * (sum(w * delta * log_t) - sum(w * cumhaz * log_t) +
                 sum(w * cumhaz_start * log_ts)) +
          g^2 * (sum(w * cumhaz * log_t^2) - sum(w * cumhaz_start * log_ts^2))

  # Design augmented with the alpha intercept column.
  x_tilde <- if (p_cov > 0) cbind(1, x) else matrix(1, nrow = n, ncol = 1L)

  hq <- crossprod(x_tilde, wA * x_tilde)             # (alpha,beta) block
  v  <- g * as.numeric(crossprod(x_tilde, wD))        # psi cross (alpha,beta)

  # Assemble p x p with order (alpha = 1, psi = 2, beta = 3..p).
  idx_q <- if (p_cov > 0) c(1L, 3:p) else 1L
  hess <- matrix(0, p, p)
  hess[idx_q, idx_q] <- hq
  hess[2L, idx_q] <- v
  hess[idx_q, 2L] <- v
  hess[2L, 2L] <- hpp
  dimnames(hess) <- list(names(theta), names(theta))
  hess
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-weibull-hessian.R")'`
Expected: PASS — all cross-check + NULL-fallback tests.

- [ ] **Step 5: Commit**

```bash
git add R/likelihood-weibull.R tests/testthat/test-weibull-hessian.R
git commit -m "feat: analytic Hessian for Weibull (internal scale, event + right-censored)"
```

---

## Task 2: Wire analytic Hessian into the Weibull optimizer

**Files:**
- Modify: `R/likelihood-weibull.R` (`.hzr_optim_weibull()` — the `.hzr_optim_generic` call ~L515)
- Test: `tests/testthat/test-weibull-hessian.R` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-weibull-hessian.R`:

```r
test_that("weibull fit vcov uses the analytic Hessian (matches numDeriv+delta)", {
  skip_if_not_installed("numDeriv")
  set.seed(25)
  n <- 500
  z <- rnorm(n)
  time <- rweibull(n, shape = 1.5, scale = exp(-0.5 * z / 1.5)) + 0.01
  status <- rbinom(n, 1, 0.85)
  fit <- hazard(
    survival::Surv(time, status) ~ z,
    data = data.frame(time, status, z),
    dist = "weibull", theta = c(mu = 1, nu = 1, z = 0), fit = TRUE
  )
  # Reference: numDeriv internal-scale Hessian -> invert -> delta-method J.
  th <- fit$fit$theta              # natural (mu, nu, beta)
  mu <- th[1]; nu <- th[2]; beta <- th[-(1:2)]
  alpha <- unname(nu * log(mu)); psi <- unname(log(nu))
  phi <- c(alpha, psi, unname(beta))
  obj <- function(p) .wb_obj_internal(p, time, status, x = cbind(z))
  v_int <- solve(numDeriv::hessian(obj, phi))
  pdim <- length(phi)
  jmat <- diag(pdim)
  jmat[1, 1] <- mu / nu
  jmat[1, 2] <- -mu * alpha / nu
  jmat[2, 2] <- nu
  v_ref <- jmat %*% v_int %*% t(jmat)
  expect_equal(unname(fit$fit$vcov), unname(v_ref), tolerance = 1e-3)
  expect_true(isTRUE(fit$fit$pd))
})

test_that("weibull SEs are invariant to covariate rescaling", {
  set.seed(26)
  n <- 600
  z <- rnorm(n)
  time <- rweibull(n, shape = 1.3, scale = exp(-0.4 * z / 1.3)) + 0.01
  status <- rbinom(n, 1, 0.9)
  f1 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z),
               dist = "weibull", theta = c(mu = 1, nu = 1, z = 0), fit = TRUE)
  f2 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z / 100),
               dist = "weibull", theta = c(mu = 1, nu = 1, z = 0), fit = TRUE)
  se1 <- sqrt(diag(f1$fit$vcov))
  se2 <- sqrt(diag(f2$fit$vcov))
  # mu and nu SEs are invariant to covariate rescaling.
  expect_equal(unname(se1[1]), unname(se2[1]), tolerance = 1e-2)
  expect_equal(unname(se1[2]), unname(se2[2]), tolerance = 1e-2)
})
```

- [ ] **Step 2: Run the tests to verify the wiring is needed**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-weibull-hessian.R")'`
Note: the vcov-match test may already pass (numDeriv is currently the source and trivially matches its own inverse+delta). Proceed to Step 3 to make the analytic Hessian the source.

- [ ] **Step 3: Pass `hessian_fn` from the Weibull optimizer**

In `R/likelihood-weibull.R`, in `.hzr_optim_weibull()`, immediately before the
`result <- .hzr_optim_generic(` call (~L515), add a closure that captures the
fit's data/weights (note: weights handled exactly as the existing `logl_internal`
does via the enclosing scope — pass the `weights` formal through):

```r
  hessian_fn <- function(par) {
    .hzr_hessian_weibull_internal(
      theta = par, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x, weights = weights
    )
  }
```

Then add `hessian_fn = hessian_fn,` as an argument to the `.hzr_optim_generic(...)`
call (e.g. immediately after `use_bounds  = FALSE`). Leave the back-transform and
the delta-method `J` block (~L527-545) unchanged — they operate on `result$vcov`
regardless of how it was produced.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-weibull-hessian.R")'`
Expected: PASS — all tests.

- [ ] **Step 5: Run the Weibull regression suites (no SE drift)**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat", filter = "weibull|gradient-weibull|weights|parity-core|wald|repeating-events")'`
Expected: 0 failures. (Analytic and numDeriv internal Hessians agree to ~1e-4, so
existing SE/vcov assertions stay within tolerance. If a tighter assertion fails,
STOP and report it — do not loosen the existing test.)

- [ ] **Step 6: Commit**

```bash
git add R/likelihood-weibull.R tests/testthat/test-weibull-hessian.R
git commit -m "feat: use analytic Hessian as Weibull SE source; add scale-invariance test"
```

---

## Task 3: NEWS entry + full check

**Files:**
- Modify: `NEWS.md`

- [ ] **Step 1: Add the NEWS entry**

In `NEWS.md`, under `# TemporalHazard 1.1.0.9000 (development version)`, in the
`## Improvements` subsection (which already has the exponential analytic-Hessian
bullet), add:

```markdown
* **Analytic Hessian for Weibull standard errors (Phase 7c, Layer 2).**
  The Weibull distribution now computes its post-fit Hessian in closed form on
  the internal `(alpha, psi, beta)` optimization scale (then mapped to the
  natural scale by the existing delta method) rather than numerically, giving
  more accurate standard errors. Covers event + right-censored data (including
  counting-process start times); left/interval-censored fits fall back to the
  numerical Hessian.
```

- [ ] **Step 2: Run the full test suite**

Run: `NOT_CRAN=true Rscript -e 'devtools::test()'`
Expected: 0 failures. (Hessian-related warnings from deliberately borderline test
fits are expected; only failures matter.)

- [ ] **Step 3: Run R CMD check**

Run: `Rscript -e 'devtools::check(args = c("--no-manual"), quiet = TRUE)'`
Expected: 0 errors / 0 warnings; only the pre-existing single NOTE if any. No new
note from `.hzr_hessian_weibull_internal`.

- [ ] **Step 4: Commit**

```bash
git add NEWS.md
git commit -m "docs: NEWS entry for Weibull analytic Hessian (Layer 2 PR-3)"
```

---

## Self-Review Notes

- **Spec coverage (Layer-2 PR-3):** analytic Weibull Hessian on the internal
  scale matching the optimizer's objective (T1), wired as prod input via the
  existing `hessian_fn` hook (T2), coverage contract = mirror gradient with
  numDeriv fallback on left/interval (T1 NULL path), analytic-vs-numDeriv cross-
  check gate incl. counting-process start (T1), end-to-end natural-scale vcov
  match + scale-invariance (T2), NEWS (T3).
- **Why internal scale:** `.hzr_optim_generic` computes/inverts the Hessian on
  the internal `phi` scale; `.hzr_optim_weibull` then applies the delta-method
  `J`. The analytic Hessian must therefore be the internal-scale objective
  Hessian, matching `numDeriv::hessian` exactly so the delta-J output is correct.
- **Foundation:** depends on PR #58/#59 (Weibull event-hazard + weighted-gradient
  fixes) being on `dev`, so `logl_internal == .hzr_logl_weibull(natural(phi))`
  for event+right and the cross-check reconstruction is valid.
- **No regression for other families:** only `.hzr_optim_weibull` gains a
  `hessian_fn`; exponential keeps its own; loglog/lognorm/multiphase still use
  numDeriv (no `hessian_fn`).
- **Type consistency:** `.hzr_hessian_weibull_internal(theta, time, status,
  time_lower, time_upper, x, weights)` matches its call site and the test calls;
  returns a p x p matrix or `NULL`.
