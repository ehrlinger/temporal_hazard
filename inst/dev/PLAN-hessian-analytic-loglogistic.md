# Hessian Stability — Layer 2 PR-4 (Log-logistic Analytic Hessian) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an analytic Hessian as the production standard-error input for the
log-logistic distribution, reusing the `hessian_fn` hook from PR-2.

**Architecture:** `.hzr_optim_loglogistic()` fits directly on
`theta = (log alpha, log beta, beta_coef)` via `.hzr_optim_generic` (no internal
reparameterization, no delta-method step — like the exponential). So the analytic
Hessian is simply the objective (negative log-likelihood) Hessian on that scale,
matching `numDeriv::hessian(objective)`. Coverage mirrors the analytic gradient:
event + right-censored only; return `NULL` on any left/interval-censored row
(`status %in% c(-1, 2)`) so the optimizer falls back to numDeriv. There is **no**
counting-process start handling in the log-logistic likelihood.

**Tech Stack:** R, testthat (edition 3), roxygen2, devtools. Base-R numerics;
cross-check tests use `numDeriv` (Suggests; guard with `skip_if_not_installed`).
Write all test code WITHOUT compound `a; b` statement lines (CI `semicolon_linter`).

**Scope:** PR-4 of the Layer-2 sequence in `inst/dev/HESSIAN-STABILITY-DESIGN.md`.
Log-logistic only. Plumbing exists (PR-2). Branch: `dev` → v1.1.0.

---

## Mathematical derivation

Parameters `theta = (a, b, beta_coef)` with `a = log(alpha)`, `b = log(beta)`,
`beta = exp(b)`. Per row, `eta_i = x_i . beta_coef`, `L_i = log(t_i)`,
`term_i = exp(a + beta * L_i + eta_i) = alpha t_i^beta exp(eta_i)`,
`p_i = term_i / (1 + term_i)`, `delta_i = (status_i == 1)`.

Log-likelihood (event + right; events weighted x2 inside the `log(1+term)` term):

    LL = sum_i w_i * delta_i * (a + b + (beta-1) L_i + eta_i)
         - sum_i w_i * (1 + delta_i) * log(1 + term_i)

Gradient building block is `d log(1+term)/d theta = p * d log(term)/d theta`,
with `d log(term)/d theta = (1, beta*L_i, x_i)` for `(a, b, beta_coef)` (note the
`b` entry is `beta*L` because `beta = exp(b)`), and `d p/d theta =
p(1-p) * d log(term)/d theta`. The only non-`log(1+term)` nonlinearity is the
`(beta-1) L` event term in `b`.

**Objective (negative log-likelihood) Hessian:** with the per-row design vector
`u_i = (1, beta*L_i, x_i)` (length p = 2 + p_cov), `m_i = w_i (1 + delta_i)`, and
`D_i = m_i p_i (1 - p_i)`:

    H_obj = U' diag(D) U                                  (U has rows u_i)
    plus an extra scalar added at the (b, b) = (2, 2) entry:
      beta * sum_i w_i L_i * ( (1 + delta_i) p_i - delta_i )

Derivation of the (2,2) extra: the `b` second derivative collects
`d^2/db^2 [ sum w delta (e^b - 1) L ] = beta * sum w delta L` (objective sign:
`- beta*sum(w*delta*L)`) and the `d^2 log(term)/db^2 = beta*L` piece inside
`log(1+term)` (`+ beta * sum(m*p*L)` in the objective). Net
`beta * sum(w*L*((1+delta)*p - delta))`. The `beta^2 L^2` curvature is already in
`U' diag(D) U` at (2,2). The numDeriv cross-check (Task 1) validates the algebra.

**Numerical guard:** use `log_t <- log(pmax(time, .Machine$double.xmin))` so a
legal right-censored `time = 0` row (where `term -> 0`, `D -> 0`) does not produce
`0 * -Inf = NaN` (same guard the Weibull Hessian uses).

---

## File Structure

- **Modify** `R/likelihood-loglogistic.R` — add `.hzr_hessian_loglogistic()`;
  wire into `.hzr_optim_loglogistic()` via a `hessian_fn` closure.
- **Test** `tests/testthat/test-loglogistic-hessian.R` (new).
- **Modify** `NEWS.md`.

---

## Task 1: `.hzr_hessian_loglogistic()` analytic Hessian

**Files:**
- Modify: `R/likelihood-loglogistic.R` (add function before `.hzr_optim_loglogistic`)
- Test: `tests/testthat/test-loglogistic-hessian.R` (create)

- [ ] **Step 1: Write the failing cross-check tests**

Create `tests/testthat/test-loglogistic-hessian.R`:

```r
test_that(".hzr_hessian_loglogistic matches numDeriv (no covariates)", {
  skip_if_not_installed("numDeriv")
  set.seed(31)
  n <- 300
  time <- rexp(n, 0.6) + 0.01
  status <- rbinom(n, 1, 0.75)
  theta <- c(log_alpha = log(0.5), log_beta = log(1.4))
  obj <- function(th) -.hzr_logl_loglogistic(th, time, status)
  h_an <- .hzr_hessian_loglogistic(theta, time, status)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_loglogistic matches numDeriv (covariates + weights)", {
  skip_if_not_installed("numDeriv")
  set.seed(32)
  n <- 400
  x <- cbind(a = rnorm(n), b = rbinom(n, 1, 0.4))
  eta <- as.numeric(x %*% c(0.3, -0.4))
  time <- rexp(n, rate = exp(-eta) * 0.5) + 0.01
  status <- rbinom(n, 1, 0.8)
  w <- runif(n, 0.5, 2)
  theta <- c(log_alpha = log(0.5), log_beta = log(1.3), 0.3, -0.4)
  obj <- function(th) -.hzr_logl_loglogistic(th, time, status, x = x, weights = w)
  h_an <- .hzr_hessian_loglogistic(theta, time, status, x = x, weights = w)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_loglogistic handles right-censored time = 0 rows", {
  skip_if_not_installed("numDeriv")
  set.seed(33)
  n <- 60
  time <- c(0, rexp(n - 1, 0.6) + 0.01)
  status <- c(0L, rbinom(n - 1, 1, 0.7))
  theta <- c(log_alpha = log(0.6), log_beta = log(1.2))
  h_an <- .hzr_hessian_loglogistic(theta, time, status)
  expect_true(all(is.finite(h_an)))
  obj <- function(th) -.hzr_logl_loglogistic(th, time, status)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_loglogistic returns NULL for left/interval censoring", {
  set.seed(34)
  n <- 40
  time <- rexp(n, 0.5) + 0.01
  status <- rep(1L, n)
  status[1:5] <- 2L
  expect_null(.hzr_hessian_loglogistic(
    c(log_alpha = 0, log_beta = 0), time, status,
    time_lower = time, time_upper = time + 0.5))
  status2 <- rep(1L, n)
  status2[1:5] <- -1L
  expect_null(.hzr_hessian_loglogistic(
    c(log_alpha = 0, log_beta = 0), time, status2, time_upper = time + 0.5))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-loglogistic-hessian.R")'`
Expected: FAIL — `could not find function ".hzr_hessian_loglogistic"`.

- [ ] **Step 3: Implement `.hzr_hessian_loglogistic()`**

In `R/likelihood-loglogistic.R`, immediately before `.hzr_optim_loglogistic()`, add:

```r
#' Analytic Hessian of the log-logistic negative log-likelihood
#'
#' Returns the Hessian of the objective (negative log-likelihood) on the internal
#' \code{theta = (log alpha, log beta, beta_coef)} scale, matching
#' \code{numDeriv::hessian(objective)} so it can be used directly by
#' \code{.hzr_optim_generic()}.  Coverage mirrors the analytic gradient: event +
#' right-censored rows only; returns \code{NULL} for any left/interval-censored
#' row (\code{status \%in\% c(-1, 2)}) to request the numerical-Hessian fallback.
#'
#' With \code{term = alpha t^beta e^eta}, \code{p = term/(1+term)},
#' \code{u = (1, beta*log t, x)} and \code{D = w (1 + delta) p (1 - p)}, the
#' objective Hessian is \code{U' diag(D) U} plus a scalar
#' \code{beta * sum(w * log t * ((1+delta) p - delta))} added at the (log beta,
#' log beta) entry.
#'
#' @noRd
.hzr_hessian_loglogistic <- function(
    theta, time, status,
    time_lower = NULL, time_upper = NULL, x = NULL, weights = NULL) {

  n <- length(time)
  if (is.null(weights)) weights <- rep(1, n)

  # Coverage contract: closed form is event + right censored only.
  if (any(status %in% c(-1, 2))) {
    return(NULL)
  }

  log_alpha <- theta[1]
  log_beta <- theta[2]
  beta <- exp(log_beta)
  p <- length(theta)
  p_cov <- p - 2L
  has_cov <- p_cov > 0L && !is.null(x)
  beta_coef <- if (has_cov) theta[3:p] else numeric(0)
  eta <- if (has_cov) as.numeric(x %*% beta_coef) else rep(0, n)

  log_t <- log(pmax(time, .Machine$double.xmin))
  term <- exp(log_alpha + beta * log_t + eta)
  pr <- term / (1 + term)
  delta <- as.numeric(status == 1)

  # Per-row design u = (1, beta*log t, x); diagonal weight D = w(1+delta)p(1-p).
  d_wt <- weights * (1 + delta) * pr * (1 - pr)
  u <- cbind(1, beta * log_t)
  if (has_cov) u <- cbind(u, x)

  hess <- crossprod(u, d_wt * u)
  hess[2L, 2L] <- hess[2L, 2L] +
    beta * sum(weights * log_t * ((1 + delta) * pr - delta))
  dimnames(hess) <- list(names(theta), names(theta))
  hess
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-loglogistic-hessian.R")'`
Expected: PASS — all cross-check + edge + NULL tests.

IF any cross-check FAILS (difference > 1e-4): do NOT tweak tolerances or
guess-patch. Report DONE_WITH_CONCERNS/BLOCKED with the analytic-vs-numDeriv
per-entry differences so the controller can correct the derivation.

- [ ] **Step 5: Commit**

```bash
git add R/likelihood-loglogistic.R tests/testthat/test-loglogistic-hessian.R
git commit -m "feat: analytic Hessian for log-logistic (event + right-censored)"
```

---

## Task 2: Wire analytic Hessian into the log-logistic optimizer

**Files:**
- Modify: `R/likelihood-loglogistic.R` (`.hzr_optim_loglogistic()`)
- Test: `tests/testthat/test-loglogistic-hessian.R` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-loglogistic-hessian.R`:

```r
test_that("loglogistic fit vcov uses the analytic Hessian (matches numDeriv)", {
  skip_if_not_installed("numDeriv")
  set.seed(35)
  n <- 500
  z <- rnorm(n)
  time <- rexp(n, rate = exp(-0.5 * z) * 0.4) + 0.01
  status <- rbinom(n, 1, 0.85)
  fit <- hazard(
    survival::Surv(time, status) ~ z,
    data = data.frame(time, status, z),
    dist = "loglogistic",
    theta = c(log_alpha = 0, log_beta = 0, z = 0), fit = TRUE
  )
  obj <- function(th) -.hzr_logl_loglogistic(th, time, status, x = cbind(z))
  v_nd <- solve(numDeriv::hessian(obj, fit$fit$theta))
  expect_equal(unname(fit$fit$vcov), unname(v_nd), tolerance = 1e-4)
  expect_true(isTRUE(fit$fit$pd))
})

test_that("loglogistic SEs are invariant to covariate rescaling", {
  set.seed(36)
  n <- 600
  z <- rnorm(n)
  time <- rexp(n, rate = exp(-0.4 * z) * 0.5) + 0.01
  status <- rbinom(n, 1, 0.9)
  f1 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z),
               dist = "loglogistic",
               theta = c(log_alpha = 0, log_beta = 0, z = 0), fit = TRUE)
  f2 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z / 100),
               dist = "loglogistic",
               theta = c(log_alpha = 0, log_beta = 0, z = 0), fit = TRUE)
  se1 <- sqrt(diag(f1$fit$vcov))
  se2 <- sqrt(diag(f2$fit$vcov))
  # log_alpha and log_beta SEs are invariant to covariate rescaling.
  expect_equal(unname(se1[1]), unname(se2[1]), tolerance = 1e-2)
  expect_equal(unname(se1[2]), unname(se2[2]), tolerance = 1e-2)
})
```

- [ ] **Step 2: Run the tests** — Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-loglogistic-hessian.R")'`
The vcov-match test may already pass (numDeriv is currently the source). Proceed to Step 3.

- [ ] **Step 3: Pass `hessian_fn` from the log-logistic optimizer**

In `R/likelihood-loglogistic.R`, replace `.hzr_optim_loglogistic()` with:

```r
#' @noRd
.hzr_optim_loglogistic <- function(
    time, status, time_lower = NULL, time_upper = NULL,
    x = NULL, theta_start, weights = NULL, control = list()) {
  hessian_fn <- function(par) {
    .hzr_hessian_loglogistic(
      theta = par, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x, weights = weights
    )
  }
  .hzr_optim_generic(
    logl_fn = .hzr_logl_loglogistic,
    gradient_fn = .hzr_gradient_loglogistic,
    time = time, status = status,
    time_lower = time_lower, time_upper = time_upper,
    x = x, theta_start = theta_start, weights = weights,
    control = control, use_bounds = FALSE,
    hessian_fn = hessian_fn
  )
}
```

(Read the current `.hzr_optim_loglogistic()` first; it calls `.hzr_optim_generic`
without `hessian_fn`. Replace the whole function with the version above.)

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-loglogistic-hessian.R")'`
Expected: PASS — all tests.

- [ ] **Step 5: Run the log-logistic regression suite (no SE drift)**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-loglogistic-dist.R")'`
Expected: 0 failures. (Analytic and numDeriv Hessians agree to ~1e-4, so existing
SE/vcov assertions stay within tolerance. If a tighter assertion fails, STOP and
report it — do not loosen the existing test.)

- [ ] **Step 6: Commit**

```bash
git add R/likelihood-loglogistic.R tests/testthat/test-loglogistic-hessian.R
git commit -m "feat: use analytic Hessian as log-logistic SE source; add scale-invariance test"
```

---

## Task 3: NEWS entry + full check

**Files:**
- Modify: `NEWS.md`

- [ ] **Step 1: Add the NEWS entry**

In `NEWS.md`, under `# TemporalHazard 1.1.0.9000 (development version)`, in the
`## Improvements` subsection (which already has the exponential and Weibull
analytic-Hessian bullets), add this bullet right after the Weibull one:

```markdown
* **Analytic Hessian for log-logistic standard errors (Phase 7c, Layer 2).**
  The log-logistic distribution now computes its post-fit Hessian in closed form
  on the internal `(log alpha, log beta, beta)` scale rather than numerically,
  giving more accurate standard errors. Covers event + right-censored data;
  left/interval-censored fits fall back to the numerical Hessian.
```

- [ ] **Step 2: Run the full test suite**

Run: `NOT_CRAN=true Rscript -e 'devtools::test()'`
Expected: 0 failures. (Borderline-fit Hessian warnings are expected; only failures
matter.)

- [ ] **Step 3: Run R CMD check**

Run: `Rscript -e 'devtools::check(args = c("--no-manual"), quiet = TRUE)'`
Expected: 0 errors / 0 warnings; only the pre-existing single NOTE if any. No new
note from `.hzr_hessian_loglogistic`.

- [ ] **Step 4: Commit**

```bash
git add NEWS.md
git commit -m "docs: NEWS entry for log-logistic analytic Hessian (Layer 2 PR-4)"
```

---

## Self-Review Notes

- **Spec coverage (Layer-2 PR-4):** analytic log-logistic Hessian matching the
  optimizer objective (T1), wired as prod input via the existing `hessian_fn`
  hook (T2), coverage contract = mirror gradient with numDeriv fallback on
  left/interval (T1 NULL path), analytic-vs-numDeriv cross-check gate incl. the
  time=0 edge (T1), end-to-end vcov match + scale-invariance (T2), NEWS (T3).
- **Direct (no delta-J):** unlike Weibull, the log-logistic optimizer fits on the
  `(log alpha, log beta, beta_coef)` scale directly and returns the generic
  result, so the analytic Hessian matches `numDeriv::hessian` of the objective
  with no transform — the end-to-end test compares `fit$fit$vcov` to the numDeriv
  inverse directly.
- **No counting-process start:** the log-logistic likelihood does not use
  `time_lower` as an epoch start (only as a left/interval bound), so the Hessian
  has no start terms — simpler than Weibull.
- **time=0 guard:** `log(pmax(time, .Machine$double.xmin))` prevents `0 * -Inf`
  on legal right-censored time=0 rows; covered by a dedicated test.
- **No regression for other families:** only `.hzr_optim_loglogistic` gains a
  `hessian_fn`; exponential/Weibull keep theirs; lognormal/multiphase still use
  numDeriv.
- **Lint:** all test code uses one statement per line (no compound `a; b`),
  per the CI `semicolon_linter`.
