# Hessian Stability — Layer 2 PR-2 (Exponential Analytic Hessian) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an analytic Hessian as the production standard-error input for the
exponential distribution, plus the reusable `hessian_fn` plumbing in the shared
optimizer that all later families (PRs 3–6) will use.

**Architecture:** `.hzr_optim_generic()` gains an optional `hessian_fn` argument;
when supplied and it returns a matrix, that matrix (the **objective / negative
log-likelihood** Hessian, same scale as `numDeriv::hessian(objective)`) is fed to
the Layer-1 `.hzr_safe_solve()` instead of the numerical Hessian. The exponential
optimizer passes a closure computing the closed-form Hessian. Coverage mirrors the
existing analytic gradient: event + right-censored rows only; if any row is
left/interval-censored (`status ∈ {-1, 2}`) the closure returns `NULL`, falling
back to numDeriv for the whole Hessian.

**Tech Stack:** R, testthat (edition 3), roxygen2, devtools. Base-R numerics; the
cross-check tests use `numDeriv` (Suggests → guarded with `skip_if_not_installed`).

**Scope:** PR-2 of the Layer-2 sequence in `inst/dev/HESSIAN-STABILITY-DESIGN.md`.
Exponential only + the shared plumbing. PRs 3–6 (Weibull, log-logistic, log-normal,
multiphase) reuse the plumbing and are out of scope here. Branch: `dev` → v1.1.0.

---

## Mathematical derivation (exponential, event + right-censored)

Parameterization (`R/likelihood-exponential.R`): `θ = [log λ, β₁, …, β_p]`,
`η_i = xᵢ·β`, cumulative hazard `H_i = λ·t_i·exp(η_i)` with `λ = exp(log λ)`.

Log-likelihood over event (`δ=1`) + right-censored (`δ=0`) rows, weights `w_i`:

    ℓ = Σ_event w_i (log λ + η_i)  −  Σ_i w_i H_i

The `log λ + η_i` event terms are **linear** in the parameters, so their second
derivatives vanish — only the `−w_i H_i` terms carry curvature. Since
`∂H_i/∂(log λ) = H_i` and `∂H_i/∂β_j = H_i·x_ij`:

Let `X̃ = [1 | X]` (n × (p+1): a leading all-ones column for `log λ`, then the
covariate columns) and `wH_i = w_i·H_i`. Then:

    Hessian(ℓ)         = − X̃ᵀ diag(wH) X̃
    Hessian(objective) = + X̃ᵀ diag(wH) X̃      (objective = −ℓ)

`.hzr_optim_generic()` minimizes the objective and inverts the **objective**
Hessian, so the analytic function returns `+X̃ᵀ diag(wH) X̃` (positive
semidefinite — Cholesky-friendly). This is exactly the derivative of the existing
analytic gradient `X̃ᵀ(w·δ − wH)`, so the two are consistent by construction.

---

## File Structure

- **Modify** `R/optimizer.R` — add `hessian_fn = NULL` to `.hzr_optim_generic()`;
  prefer it over numDeriv when it returns a matrix.
- **Modify** `R/likelihood-exponential.R` — add `.hzr_hessian_exponential()`;
  wire it into `.hzr_optim_exponential()`.
- **Test** `tests/testthat/test-hessian-invert.R` — plumbing tests for the generic.
- **Test** `tests/testthat/test-exponential-hessian.R` (new) — analytic-vs-numDeriv
  cross-check, NULL fallback, end-to-end fit vcov, scale-invariance.
- **Modify** `NEWS.md` — entry under `# TemporalHazard 1.1.0.9000`.

---

## Task 1: `hessian_fn` plumbing in `.hzr_optim_generic()`

**Files:**
- Modify: `R/optimizer.R` (signature ~L42–54; Hessian block ~L129–146)
- Test: `tests/testthat/test-hessian-invert.R` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-hessian-invert.R`:

```r
test_that(".hzr_optim_generic uses a supplied analytic hessian_fn", {
  set.seed(1)
  time <- rexp(100, rate = 0.5); status <- rep(1L, 100)
  marker <- matrix(42, 1, 1)
  res <- .hzr_optim_generic(
    logl_fn = .hzr_logl_exponential, gradient_fn = .hzr_gradient_exponential,
    time = time, status = status, x = NULL, theta_start = c(log_rate = 0),
    control = list(maxit = 200, reltol = 1e-8, abstol = 1e-8),
    hessian_fn = function(par) marker
  )
  expect_equal(unname(res$hessian), unname(marker))
  expect_equal(unname(res$vcov), solve(marker), tolerance = 1e-10)
})

test_that(".hzr_optim_generic falls back to numDeriv when hessian_fn returns NULL", {
  skip_if_not_installed("numDeriv")
  set.seed(1)
  time <- rexp(100, rate = 0.5); status <- rep(1L, 100)
  res <- .hzr_optim_generic(
    logl_fn = .hzr_logl_exponential, gradient_fn = .hzr_gradient_exponential,
    time = time, status = status, x = NULL, theta_start = c(log_rate = 0),
    control = list(maxit = 200, reltol = 1e-8, abstol = 1e-8),
    hessian_fn = function(par) NULL
  )
  expect_true(is.matrix(res$hessian) && is.finite(res$hessian[1, 1]))
  expect_true(isTRUE(res$pd))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hessian-invert.R")'`
Expected: FAIL — `unused argument (hessian_fn = ...)`.

- [ ] **Step 3: Add the `hessian_fn` parameter**

In `R/optimizer.R`, change the signature (currently ending
`use_bounds = FALSE, lower_bounds = NULL)`) to add `hessian_fn = NULL`:

```r
.hzr_optim_generic <- function(
    logl_fn,
    gradient_fn,
    time,
    status,
    time_lower = NULL,
    time_upper = NULL,
    x = NULL,
    theta_start,
    weights = NULL,
    control = list(),
    use_bounds = FALSE,
    lower_bounds = NULL,
    hessian_fn = NULL) {
```

Also add a roxygen `@param` line (above the function, with the other `@param`
entries):

```r
#' @param hessian_fn Optional function(theta) returning the Hessian of the
#'   negative log-likelihood at \code{theta} (same scale as
#'   \code{numDeriv::hessian(objective)}).  When it returns a matrix, that matrix
#'   is used for standard errors; when \code{NULL} (the default) or when the
#'   function returns \code{NULL}, a numerical Hessian is used instead.
```

- [ ] **Step 4: Use `hessian_fn` in the post-fit Hessian block**

In `R/optimizer.R`, replace the post-fit Hessian block (from the
`# Post-fit Hessian for standard errors` comment through the
`error = function(e) NA\n  )` that closes the numDeriv `tryCatch`) with:

```r
  # Post-fit Hessian for standard errors.  Prefer the caller's analytic Hessian
  # (on the objective / negative-log-likelihood scale) when supplied; otherwise,
  # or when it declines by returning NULL, fall back to a numerical Hessian.
  hess_result <- NULL
  if (!is.null(hessian_fn)) {
    hess_result <- tryCatch(hessian_fn(result$par), error = function(e) NULL)
  }
  if (is.null(hess_result)) {
    hess_result <- tryCatch(
      {
        if (requireNamespace("numDeriv", quietly = TRUE)) {
          numDeriv::hessian(objective, result$par)
        } else {
          NA
        }
      },
      error = function(e) NA
    )
  }
```

Leave the existing `inv <- if (is.matrix(hess_result)) .hzr_safe_solve(...) ...`
block and the final `list(...)` return unchanged.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hessian-invert.R")'`
Expected: PASS (all prior tests plus the two new ones).

- [ ] **Step 6: Commit**

```bash
git add R/optimizer.R tests/testthat/test-hessian-invert.R
git commit -m "feat: add hessian_fn plumbing to .hzr_optim_generic (Layer 2)"
```

---

## Task 2: `.hzr_hessian_exponential()` analytic Hessian

**Files:**
- Modify: `R/likelihood-exponential.R` (add function near `.hzr_gradient_exponential`)
- Test: `tests/testthat/test-exponential-hessian.R` (create)

- [ ] **Step 1: Write the failing cross-check tests**

Create `tests/testthat/test-exponential-hessian.R`:

```r
test_that(".hzr_hessian_exponential matches numDeriv (no covariates)", {
  skip_if_not_installed("numDeriv")
  set.seed(2)
  n <- 300
  time <- rexp(n, rate = 0.7); status <- rbinom(n, 1, 0.7)
  theta <- c(log_rate = log(0.7))
  obj <- function(par) -.hzr_logl_exponential(par, time, status, x = NULL)
  H_an <- .hzr_hessian_exponential(theta, time, status, x = NULL)
  H_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_exponential matches numDeriv (covariates + weights)", {
  skip_if_not_installed("numDeriv")
  set.seed(3)
  n <- 400
  x <- cbind(a = rnorm(n), b = rbinom(n, 1, 0.4))
  eta <- as.numeric(x %*% c(0.3, -0.5))
  time <- rexp(n, rate = exp(log(0.5) + eta)); status <- rbinom(n, 1, 0.8)
  w <- runif(n, 0.5, 2)
  theta <- c(log_rate = log(0.5), 0.3, -0.5)
  obj <- function(par) -.hzr_logl_exponential(par, time, status, x = x, weights = w)
  H_an <- .hzr_hessian_exponential(theta, time, status, x = x, weights = w)
  H_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_exponential returns NULL for left/interval censoring", {
  set.seed(4)
  n <- 50
  time <- rexp(n, 0.5); status <- rep(1L, n); status[1:5] <- 2L
  tl <- time; tu <- time + 0.5
  expect_null(.hzr_hessian_exponential(c(log_rate = 0), time, status,
                                       time_lower = tl, time_upper = tu, x = NULL))
  status2 <- rep(1L, n); status2[1:5] <- -1L
  expect_null(.hzr_hessian_exponential(c(log_rate = 0), time, status2,
                                       time_upper = tu, x = NULL))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-exponential-hessian.R")'`
Expected: FAIL — `could not find function ".hzr_hessian_exponential"`.

- [ ] **Step 3: Implement `.hzr_hessian_exponential()`**

In `R/likelihood-exponential.R`, immediately after `.hzr_gradient_exponential()`
(before `.hzr_optim_exponential()`), add:

```r
#' Analytic Hessian of the exponential negative log-likelihood
#'
#' Returns the Hessian of the **objective** (negative log-likelihood) on the
#' internal parameter scale, matching \code{numDeriv::hessian(objective)} so it
#' can be used directly by \code{.hzr_optim_generic()}.  Coverage mirrors the
#' analytic gradient: closed form for event + right-censored rows only.  When any
#' row is left- or interval-censored (\code{status \%in\% c(-1, 2)}) it returns
#' \code{NULL} to request the numerical-Hessian fallback.
#'
#' Derivation: with \eqn{H_i = \lambda t_i e^{\eta_i}} and design matrix
#' \eqn{\tilde X = [1 | X]}, the log-likelihood Hessian is
#' \eqn{-\tilde X^\top \mathrm{diag}(w H) \tilde X}; the objective Hessian is its
#' negation, \eqn{+\tilde X^\top \mathrm{diag}(w H) \tilde X}.
#'
#' @noRd
.hzr_hessian_exponential <- function(
    theta, time, status,
    time_lower = NULL, time_upper = NULL, x = NULL, weights = NULL) {

  n <- length(time)
  if (is.null(weights)) weights <- rep(1, n)

  # Mirror the analytic gradient's coverage: closed form is event + right
  # censored only.  Decline (NULL -> numDeriv fallback) on left/interval rows.
  if (any(status %in% c(-1, 2))) {
    return(NULL)
  }

  log_lambda <- theta[1]
  lambda <- exp(log_lambda)
  if (!is.null(x)) {
    beta <- theta[2:length(theta)]
    eta <- as.numeric(x %*% beta)
  } else {
    eta <- rep(0, n)
  }

  # Weighted cumulative hazard per row.
  wH <- weights * lambda * time * exp(eta)

  # Design matrix augmented with the log(lambda) "intercept" column.
  x_tilde <- if (is.null(x)) matrix(1, nrow = n, ncol = 1L) else cbind(1, x)

  # Objective (negative log-likelihood) Hessian: + X~' diag(wH) X~.
  hess <- crossprod(x_tilde, wH * x_tilde)
  dimnames(hess) <- list(names(theta), names(theta))
  hess
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-exponential-hessian.R")'`
Expected: PASS — all cross-check and NULL-fallback tests pass.

- [ ] **Step 5: Commit**

```bash
git add R/likelihood-exponential.R tests/testthat/test-exponential-hessian.R
git commit -m "feat: analytic Hessian for exponential (event + right-censored)"
```

---

## Task 3: Wire analytic Hessian into the exponential optimizer

**Files:**
- Modify: `R/likelihood-exponential.R` (`.hzr_optim_exponential()`)
- Test: `tests/testthat/test-exponential-hessian.R` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-exponential-hessian.R`:

```r
test_that("exponential fit vcov uses the analytic Hessian (matches numDeriv)", {
  skip_if_not_installed("numDeriv")
  set.seed(5)
  n <- 500
  z <- rnorm(n)
  time <- rexp(n, rate = exp(log(0.4) + 0.6 * z)); status <- rbinom(n, 1, 0.85)
  fit <- hazard(
    survival::Surv(time, status) ~ z,
    data = data.frame(time, status, z),
    dist = "exponential", theta = c(log_rate = 0, z = 0), fit = TRUE
  )
  obj <- function(par) -.hzr_logl_exponential(par, time, status, x = cbind(z))
  v_nd <- solve(numDeriv::hessian(obj, fit$fit$theta))
  expect_equal(unname(fit$fit$vcov), unname(v_nd), tolerance = 1e-4)
  expect_true(isTRUE(fit$fit$pd))
})

test_that("exponential SEs are invariant to covariate rescaling", {
  set.seed(6)
  n <- 600
  z <- rnorm(n, mean = 50, sd = 10)
  time <- rexp(n, rate = exp(log(0.3) + 0.02 * (z - 50)))
  status <- rbinom(n, 1, 0.9)
  f1 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z),
               dist = "exponential", theta = c(log_rate = 0, z = 0), fit = TRUE)
  f2 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z / 100),
               dist = "exponential", theta = c(log_rate = 0, z = 0), fit = TRUE)
  se1 <- sqrt(diag(f1$fit$vcov)); se2 <- sqrt(diag(f2$fit$vcov))
  # log_rate SE is invariant to covariate rescaling.
  expect_equal(unname(se1[1]), unname(se2[1]), tolerance = 1e-2)
  # The beta z-statistic is invariant under x -> x / 100.
  z1 <- unname(f1$fit$theta[2] / se1[2])
  z2 <- unname(f2$fit$theta[2] / se2[2])
  expect_equal(z1, z2, tolerance = 1e-2)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-exponential-hessian.R")'`
Expected: FAIL on the vcov-matches-numDeriv test only if the analytic path is not
yet wired (it currently uses numDeriv, which would also match — so this test may
already pass; the wiring test below is the definitive one). If both pass before
wiring, proceed to Step 3 anyway: wiring makes the analytic Hessian the *source*.

(Determinism note: the scale-invariance test sets its own seed; both fits use
the package's multi-start default which is irrelevant for single-distribution
fits — exponential uses a single BFGS run.)

- [ ] **Step 3: Pass `hessian_fn` from the exponential optimizer**

In `R/likelihood-exponential.R`, replace `.hzr_optim_exponential()` with:

```r
#' @noRd
.hzr_optim_exponential <- function(
    time, status, time_lower = NULL, time_upper = NULL,
    x = NULL, theta_start, weights = NULL, control = list()) {
  hessian_fn <- function(par) {
    .hzr_hessian_exponential(
      theta = par, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x, weights = weights
    )
  }
  .hzr_optim_generic(
    logl_fn = .hzr_logl_exponential,
    gradient_fn = .hzr_gradient_exponential,
    time = time, status = status,
    time_lower = time_lower, time_upper = time_upper,
    x = x, theta_start = theta_start, weights = weights,
    control = control, use_bounds = FALSE,
    hessian_fn = hessian_fn
  )
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-exponential-hessian.R")'`
Expected: PASS — all tests.

- [ ] **Step 5: Run the exponential regression suite (no SE drift)**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-exponential-dist.R")'`
Expected: 0 failures. (Analytic and numDeriv Hessians agree to ~1e-4, so any
existing SE/vcov assertions remain within tolerance. If a tighter assertion
fails, STOP and report it — do not loosen the existing test.)

- [ ] **Step 6: Commit**

```bash
git add R/likelihood-exponential.R tests/testthat/test-exponential-hessian.R
git commit -m "feat: use analytic Hessian as exponential SE source; add scale-invariance test"
```

---

## Task 4: NEWS entry + full check

**Files:**
- Modify: `NEWS.md`

- [ ] **Step 1: Add the NEWS entry**

In `NEWS.md`, under `# TemporalHazard 1.1.0.9000 (development version)`, in the
`## Improvements` subsection (the one added by the Layer-1 PR), add:

```markdown
* **Analytic Hessian for exponential standard errors (Phase 7c, Layer 2).**
  The exponential distribution now computes its post-fit Hessian in closed form
  (`X~' diag(wH) X~` over event + right-censored rows) rather than numerically,
  giving faster and more accurate standard errors. The shared optimizer gained a
  `hessian_fn` hook that analytic Hessians for the remaining families will reuse;
  left/interval-censored exponential fits fall back to the numerical Hessian.
```

- [ ] **Step 2: Run the full test suite**

Run: `NOT_CRAN=true Rscript -e 'devtools::test()'`
Expected: 0 failures across the whole suite.

- [ ] **Step 3: Run R CMD check**

Run: `Rscript -e 'devtools::check(args = c("--no-manual"), quiet = TRUE)'`
Expected: 0 errors / 0 warnings; only the pre-existing single NOTE (if any).
Confirm no new note from the new internal `.hzr_hessian_exponential` symbol.

- [ ] **Step 4: Commit**

```bash
git add NEWS.md
git commit -m "docs: NEWS entry for exponential analytic Hessian (Layer 2 PR-2)"
```

---

## Self-Review Notes

- **Spec coverage (Layer 2 / PR-2 portion of HESSIAN-STABILITY-DESIGN.md):**
  `hessian_fn` plumbing (T1), exponential analytic Hessian as prod input (T2/T3),
  coverage contract = mirror gradient with numDeriv fallback on left/interval
  (T2 NULL path), analytic-vs-numDeriv cross-check gate (T2), scale-invariance
  test (T3), NEWS (T4). PRs 3–6 (other families) intentionally absent.
- **Objective-scale consistency:** both `hessian_fn` and `numDeriv::hessian` return
  the Hessian of the *objective* (negative log-likelihood); `.hzr_safe_solve`
  inverts that to vcov. The analytic formula returns `+X̃ᵀ diag(wH) X̃` to match.
- **No regression for other families:** Weibull/log-logistic/log-normal/multiphase
  do not pass `hessian_fn`, so `hess_result` stays `NULL` → numDeriv path →
  byte-identical behavior to before. Only exponential changes its SE source.
- **Type consistency:** `.hzr_hessian_exponential(theta, time, status, time_lower,
  time_upper, x, weights)` signature matches both its call site in
  `.hzr_optim_exponential` and the test call sites; returns a matrix or `NULL`.
