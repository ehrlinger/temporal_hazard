# Hessian Stability — Layer 2 PR-5 (Log-normal Analytic Hessian) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an analytic Hessian as the production standard-error input for the
log-normal (AFT) distribution, reusing the `hessian_fn` hook from PR-2.

**Architecture:** `.hzr_optim_lognormal()` fits directly on
`theta = (mu, log_sigma, beta_coef)` via `.hzr_optim_generic` (no reparam, no
delta-method — like exponential/log-logistic), where `eta = mu + x . beta_coef`,
`sigma = exp(log_sigma)`, `z = (log t - eta)/sigma`. So the analytic Hessian is
the objective (negative log-likelihood) Hessian on that scale, matching
`numDeriv::hessian(objective)`. Coverage mirrors the analytic gradient: event +
right-censored only; return `NULL` on any left/interval-censored row
(`status %in% c(-1, 2)`). No counting-process start.

**Tech Stack:** R, testthat (edition 3), roxygen2, devtools. Base-R numerics;
cross-check tests use `numDeriv` (Suggests; `skip_if_not_installed`). All test
code one-statement-per-line (CI `semicolon_linter`).

**Scope:** PR-5 of the Layer-2 sequence in `inst/dev/HESSIAN-STABILITY-DESIGN.md`.
Log-normal only. Plumbing exists (PR-2). Branch: `dev` → v1.1.0.

---

## Mathematical derivation

Per row: `eta_i = mu + x_i . beta_coef`, `sigma = exp(log_sigma)`,
`z_i = (log t_i - eta_i)/sigma`, `delta_i = (status_i == 1)`, weight `w_i`.
Inverse Mills ratio `m_i = phi(z_i)/Phi(-z_i)` with the key identity
`dm/dz = m(m - z)`.

Log-likelihood (event + right):
- event: `-log_sigma - log t - 0.5*log(2pi) - 0.5 z^2`
- right: `log Phi(-z)`

Gradient (matches `.hzr_gradient_lognormal`): `dl/d eta = (delta*z + (1-delta)*m)/sigma`,
`dl/d log_sigma = delta*(z^2-1) + (1-delta)*m*z`. Because `eta = mu + x.beta`,
`d eta/d(mu, beta) = (1, x)` (mu is the location intercept).

**Objective (negative log-likelihood) Hessian** — per-row coefficients (negation
of the LL second derivatives), then mapped through the design `Xtilde = [1 | X]`
for the `(mu, beta)` block and the separate `log_sigma` (`s`) parameter:

    A_i = ( delta + (1-delta) * m*(m - z) ) / sigma^2                 # (eta, eta)
    B_i = ( delta*2z + (1-delta) * m*(z*(m - z) + 1) ) / sigma        # (eta, s)
    C_i =   delta*2z^2 + (1-delta) * m*z*(z*(m - z) + 1)              # (s, s)

    H_obj[(mu,beta),(mu,beta)] = Xtilde' diag(w*A) Xtilde
    H_obj[s, (mu,beta)]        = Xtilde' (w*B)
    H_obj[s, s]                = sum(w*C)

Assemble into p x p (p = 2 + p_cov) at order (mu=1, s=2, beta=3..p): the Xtilde
block at `idx_q = c(1, 3:p)`, `s` at index 2.

(Derivation uses `dz/d eta = -1/sigma`, `dz/ds = -z`, `dm/dz = m(m-z)`. The event
terms come from `-0.5 z^2 - s`; the right-censored terms from `log Phi(-z)`. The
numDeriv cross-check (Task 1) is the gate that validates the algebra.)

**Mills clamp:** compute `m` exactly as the gradient does — on the log scale with
`m <- pmin(exp(dnorm(z, log=TRUE) - pnorm(-z, log.p=TRUE)), 1e6)` — for numerical
stability and consistency with the score. (At reasonable `theta` the clamp never
binds, so analytic == numDeriv on the cross-check data.)

**No time=0 guard needed:** the log-normal LL rejects `status %in% c(1, 0) &
time <= 0` (returns `Inf`), so right-censored `time = 0` is infeasible and never
reaches the Hessian — unlike Weibull/log-logistic. (Confirmed in the family
gradient audit.)

---

## File Structure

- **Modify** `R/likelihood-lognormal.R` — add `.hzr_hessian_lognormal()`; wire
  into `.hzr_optim_lognormal()` via a `hessian_fn` closure.
- **Test** `tests/testthat/test-lognormal-hessian.R` (new).
- **Modify** `NEWS.md`.

---

## Task 1: `.hzr_hessian_lognormal()` analytic Hessian

**Files:**
- Modify: `R/likelihood-lognormal.R` (add function before `.hzr_optim_lognormal`)
- Test: `tests/testthat/test-lognormal-hessian.R` (create)

- [ ] **Step 1: Write the failing cross-check tests**

Create `tests/testthat/test-lognormal-hessian.R`:

```r
test_that(".hzr_hessian_lognormal matches numDeriv (no covariates)", {
  skip_if_not_installed("numDeriv")
  set.seed(51)
  n <- 300
  time <- exp(rnorm(n, 0.2, 0.8))
  status <- rbinom(n, 1, 0.75)
  theta <- c(mu = 0.2, log_sigma = log(0.8))
  obj <- function(th) -.hzr_logl_lognormal(th, time, status)
  h_an <- .hzr_hessian_lognormal(theta, time, status)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_lognormal matches numDeriv (covariates + weights)", {
  skip_if_not_installed("numDeriv")
  set.seed(52)
  n <- 400
  x <- cbind(a = rnorm(n), b = rbinom(n, 1, 0.4))
  eta <- 0.2 + as.numeric(x %*% c(0.3, -0.4))
  time <- exp(rnorm(n, eta, 0.7))
  status <- rbinom(n, 1, 0.8)
  w <- runif(n, 0.5, 2)
  theta <- c(mu = 0.2, log_sigma = log(0.7), 0.3, -0.4)
  obj <- function(th) -.hzr_logl_lognormal(th, time, status, x = x, weights = w)
  h_an <- .hzr_hessian_lognormal(theta, time, status, x = x, weights = w)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_lognormal matches numDeriv (heavy right-censoring)", {
  skip_if_not_installed("numDeriv")
  set.seed(53)
  n <- 300
  time <- exp(rnorm(n, 0.0, 0.6))
  status <- rbinom(n, 1, 0.3)
  theta <- c(mu = 0.1, log_sigma = log(0.6))
  obj <- function(th) -.hzr_logl_lognormal(th, time, status)
  h_an <- .hzr_hessian_lognormal(theta, time, status)
  h_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(h_an), unname(h_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_lognormal returns NULL for left/interval censoring", {
  set.seed(54)
  n <- 40
  time <- exp(rnorm(n, 0, 0.5))
  status <- rep(1L, n)
  status[1:5] <- 2L
  expect_null(.hzr_hessian_lognormal(
    c(mu = 0, log_sigma = 0), time, status,
    time_lower = time, time_upper = time + 0.5))
  status2 <- rep(1L, n)
  status2[1:5] <- -1L
  expect_null(.hzr_hessian_lognormal(
    c(mu = 0, log_sigma = 0), time, status2, time_upper = time + 0.5))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-lognormal-hessian.R")'`
Expected: FAIL — `could not find function ".hzr_hessian_lognormal"`.

- [ ] **Step 3: Implement `.hzr_hessian_lognormal()`**

In `R/likelihood-lognormal.R`, immediately before `.hzr_optim_lognormal()`, add:

```r
#' Analytic Hessian of the log-normal negative log-likelihood
#'
#' Returns the Hessian of the objective (negative log-likelihood) on the internal
#' \code{theta = (mu, log_sigma, beta_coef)} scale, matching
#' \code{numDeriv::hessian(objective)} so it can be used directly by
#' \code{.hzr_optim_generic()}.  Coverage mirrors the analytic gradient: event +
#' right-censored rows only; returns \code{NULL} for any left/interval-censored
#' row (\code{status \%in\% c(-1, 2)}) to request the numerical-Hessian fallback.
#'
#' With \code{z = (log t - eta)/sigma}, \code{eta = mu + x beta_coef}, and the
#' inverse Mills ratio \code{m = phi(z)/Phi(-z)} (\code{dm/dz = m(m - z)}), the
#' objective Hessian assembles per-row coefficients A (eta,eta), B (eta,log_sigma),
#' C (log_sigma,log_sigma) through the design \code{Xtilde = [1 | X]} for the
#' \code{(mu, beta_coef)} block.
#'
#' @noRd
.hzr_hessian_lognormal <- function(
    theta, time, status,
    time_lower = NULL, time_upper = NULL, x = NULL, weights = NULL) {

  n <- length(time)
  if (is.null(weights)) weights <- rep(1, n)

  # Coverage contract: closed form is event + right censored only.
  if (any(status %in% c(-1, 2))) {
    return(NULL)
  }

  mu <- theta[1]
  log_sigma <- theta[2]
  sigma <- exp(log_sigma)
  p <- length(theta)
  p_cov <- p - 2L
  has_cov <- p_cov > 0L && !is.null(x)
  beta_coef <- if (has_cov) theta[3:p] else numeric(0)
  eta <- if (has_cov) mu + as.numeric(x %*% beta_coef) else rep(mu, n)

  z <- (log(time) - eta) / sigma
  # Inverse Mills ratio on the log scale, clamped like the analytic gradient.
  m <- pmin(exp(dnorm(z, log = TRUE) - pnorm(-z, log.p = TRUE)), 1e6)

  delta <- as.numeric(status == 1)
  cens <- 1 - delta

  # Per-row objective Hessian coefficients (negation of the LL second derivs).
  a_coef <- (delta + cens * m * (m - z)) / sigma^2                 # (eta, eta)
  b_coef <- (delta * 2 * z + cens * m * (z * (m - z) + 1)) / sigma # (eta, log_sigma)
  c_coef <- delta * 2 * z^2 + cens * m * z * (z * (m - z) + 1)     # (log_sigma)^2

  x_tilde <- if (has_cov) cbind(1, x) else matrix(1, nrow = n, ncol = 1L)

  hqq <- crossprod(x_tilde, (weights * a_coef) * x_tilde)
  v <- as.numeric(crossprod(x_tilde, weights * b_coef))
  hss <- sum(weights * c_coef)

  # Assemble p x p with order (mu = 1, log_sigma = 2, beta_coef = 3..p).
  idx_q <- if (has_cov) c(1L, 3:p) else 1L
  hess <- matrix(0, p, p)
  hess[idx_q, idx_q] <- hqq
  hess[2L, idx_q] <- v
  hess[idx_q, 2L] <- v
  hess[2L, 2L] <- hss
  dimnames(hess) <- list(names(theta), names(theta))
  hess
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-lognormal-hessian.R")'`
Expected: PASS — all cross-check + NULL tests.

IF any cross-check FAILS (> 1e-4): do NOT tweak tolerances or guess-patch.
Report DONE_WITH_CONCERNS/BLOCKED with the analytic-vs-numDeriv per-entry
differences so the controller can correct the derivation (the Mills-ratio second
derivatives are the highest-risk part).

- [ ] **Step 5: Commit**

```bash
git add R/likelihood-lognormal.R tests/testthat/test-lognormal-hessian.R
git commit -m "feat: analytic Hessian for log-normal (event + right-censored)"
```

---

## Task 2: Wire analytic Hessian into the log-normal optimizer

**Files:**
- Modify: `R/likelihood-lognormal.R` (`.hzr_optim_lognormal()`)
- Test: `tests/testthat/test-lognormal-hessian.R` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-lognormal-hessian.R`:

```r
test_that("lognormal fit vcov uses the analytic Hessian (matches numDeriv)", {
  skip_if_not_installed("numDeriv")
  set.seed(55)
  n <- 500
  z <- rnorm(n)
  time <- exp(rnorm(n, 0.2 + 0.5 * z, 0.7))
  status <- rbinom(n, 1, 0.85)
  fit <- hazard(
    survival::Surv(time, status) ~ z,
    data = data.frame(time, status, z),
    dist = "lognormal",
    theta = c(mu = 0, log_sigma = 0, z = 0), fit = TRUE
  )
  obj <- function(th) -.hzr_logl_lognormal(th, time, status, x = cbind(z))
  v_nd <- solve(numDeriv::hessian(obj, fit$fit$theta))
  expect_equal(unname(fit$fit$vcov), unname(v_nd), tolerance = 1e-4)
  expect_true(isTRUE(fit$fit$pd))
})

test_that("lognormal SEs are invariant to covariate rescaling", {
  set.seed(56)
  n <- 600
  z <- rnorm(n)
  time <- exp(rnorm(n, 0.1 + 0.4 * z, 0.6))
  status <- rbinom(n, 1, 0.9)
  f1 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z),
               dist = "lognormal",
               theta = c(mu = 0, log_sigma = 0, z = 0), fit = TRUE)
  f2 <- hazard(survival::Surv(time, status) ~ z,
               data = data.frame(time, status, z = z / 100),
               dist = "lognormal",
               theta = c(mu = 0, log_sigma = 0, z = 0), fit = TRUE)
  se1 <- sqrt(diag(f1$fit$vcov))
  se2 <- sqrt(diag(f2$fit$vcov))
  # mu and log_sigma SEs are invariant to covariate rescaling.
  expect_equal(unname(se1[1]), unname(se2[1]), tolerance = 1e-2)
  expect_equal(unname(se1[2]), unname(se2[2]), tolerance = 1e-2)
})
```

- [ ] **Step 2: Run the tests** — Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-lognormal-hessian.R")'`
The vcov-match test may already pass (numDeriv is currently the source). Proceed to Step 3.

- [ ] **Step 3: Pass `hessian_fn` from the optimizer** — In `R/likelihood-lognormal.R`, replace `.hzr_optim_lognormal()` with:

```r
#' @noRd
.hzr_optim_lognormal <- function(
    time, status, time_lower = NULL, time_upper = NULL,
    x = NULL, theta_start, weights = NULL, control = list()) {
  hessian_fn <- function(par) {
    .hzr_hessian_lognormal(
      theta = par, time = time, status = status,
      time_lower = time_lower, time_upper = time_upper,
      x = x, weights = weights
    )
  }
  .hzr_optim_generic(
    logl_fn = .hzr_logl_lognormal,
    gradient_fn = .hzr_gradient_lognormal,
    time = time, status = status,
    time_lower = time_lower, time_upper = time_upper,
    x = x, theta_start = theta_start, weights = weights,
    control = control, use_bounds = FALSE,
    hessian_fn = hessian_fn
  )
}
```

(Read the current `.hzr_optim_lognormal()` first; it calls `.hzr_optim_generic`
without `hessian_fn`. Replace the whole function with the version above.)

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-lognormal-hessian.R")'`
Expected: PASS — all tests.

- [ ] **Step 5: Run the log-normal regression suite (no SE drift)**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-lognormal-dist.R")'`
Expected: 0 failures. (If a tighter SE assertion fails, STOP and report — do not
loosen the existing test.)

- [ ] **Step 6: Commit**

```bash
git add R/likelihood-lognormal.R tests/testthat/test-lognormal-hessian.R
git commit -m "feat: use analytic Hessian as log-normal SE source; add scale-invariance test"
```

---

## Task 3: NEWS entry + full check

**Files:**
- Modify: `NEWS.md`

- [ ] **Step 1: Add the NEWS entry**

In `NEWS.md`, under `# TemporalHazard 1.1.0.9000 (development version)`, in the
`## Improvements` subsection (after the log-logistic analytic-Hessian bullet), add:

```markdown
* **Analytic Hessian for log-normal standard errors (Phase 7c, Layer 2).**
  The log-normal distribution now computes its post-fit Hessian in closed form
  on the internal `(mu, log_sigma, beta_coef)` scale rather than numerically,
  giving more accurate standard errors. Covers event + right-censored data;
  left/interval-censored fits fall back to the numerical Hessian.
```

- [ ] **Step 2: Run the full test suite**

Run: `NOT_CRAN=true Rscript -e 'devtools::test()'`
Expected: 0 failures.

- [ ] **Step 3: Run R CMD check**

Run: `Rscript -e 'devtools::check(args = c("--no-manual"), quiet = TRUE)'`
Expected: 0 errors / 0 warnings; only the pre-existing single NOTE if any. No new
note from `.hzr_hessian_lognormal`.

- [ ] **Step 4: Commit**

```bash
git add NEWS.md
git commit -m "docs: NEWS entry for log-normal analytic Hessian (Layer 2 PR-5)"
```

---

## Self-Review Notes

- **Spec coverage (Layer-2 PR-5):** analytic log-normal Hessian matching the
  optimizer objective (T1), wired as prod input (T2), coverage contract = mirror
  gradient with numDeriv fallback on left/interval (T1 NULL path), cross-check
  gate incl. heavy right-censoring (T1), end-to-end vcov match + scale-invariance
  (T2), NEWS (T3).
- **Direct (no delta-J):** log-normal fits on `(mu, log_sigma, beta_coef)` and
  returns the generic result, so the analytic Hessian matches `numDeriv::hessian`
  of the objective directly — end-to-end test compares `fit$fit$vcov` to the
  numDeriv inverse with no transform.
- **No time=0 guard:** the LL rejects `status %in% c(1, 0) & time <= 0`, so the
  Hessian never sees `time = 0` (unlike Weibull/log-logistic).
- **Mills clamp** matches the gradient (`pmin(..., 1e6)`); never binds at the
  cross-check thetas, so analytic == numDeriv.
- **No regression for other families:** only `.hzr_optim_lognormal` gains a
  `hessian_fn`; exp/weibull/loglogistic keep theirs; multiphase still numDeriv.
- **Lint:** test code one statement per line (no compound `a; b`).
