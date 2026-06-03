# Hessian Stability — Layer 1 (PR-1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the bare `solve()` Hessian inversion with a hardened,
diagnostic-emitting helper so no fit silently returns garbage standard
errors, and surface `rcond`/`pd` diagnostics on the fit and in `summary()`.

**Architecture:** Add `.hzr_safe_solve()` (symmetrize → rcond check →
Cholesky-with-`solve()`-fallback → negative-variance guard, with named
warnings) in a new focused file `R/hessian-invert.R`. Call it from the single
shared chokepoint `.hzr_optim_generic()` in `R/optimizer.R`; thread the
returned `rcond`/`pd` into `fit_state` in `R/hazard_api.R`; print a flagged
note in `print.summary.hazard()`.

**Tech Stack:** R, testthat (edition 3), roxygen2, devtools. Base-R numerics
only (`rcond`, `chol`, `chol2inv`, `solve`); no new dependencies.

**Scope:** This is PR-1 of the 6-PR sequence in
`inst/dev/HESSIAN-STABILITY-DESIGN.md`. It covers **Layer 1 only** (universal
inversion hardening + diagnostics). The analytic Hessians (PRs 2–6) are out
of scope here and get their own plans. Branch: `dev` → v1.1.0.

---

## File Structure

- **Create** `R/hessian-invert.R` — `.hzr_safe_solve()` + the shared
  `.hzr_rcond_tol` threshold constant. Single responsibility: turn a raw
  Hessian into a vcov + diagnostics.
- **Create** `tests/testthat/test-hessian-invert.R` — unit tests for
  `.hzr_safe_solve()` against hand-built matrices.
- **Modify** `R/optimizer.R` (~L129–162) — call `.hzr_safe_solve()` instead
  of `solve()`; add `rcond`/`pd` to the returned list.
- **Modify** `R/hazard_api.R` (~L444 and ~L474) — copy `rcond`/`pd` into
  `fit_state` in both fit branches.
- **Modify** `R/hazard_api.R` (`summary.hazard` ~L1121, `print.summary.hazard`
  ~L1180) — store and print the flagged diagnostic note.
- **Create** `tests/testthat/test-hessian-stability.R` — the 13-parameter
  multiphase anchor regression.
- **Modify** `NEWS.md` — Layer-1 entry under `# TemporalHazard 1.1.0.9000`.

**Design note (surfacing policy):** Layer 1 warnings for the *new* degenerate
paths (ill-conditioned, non-positive variance, non-finite) are NOT added to
the muffle list in `.hzr_run_fit_safely()` (`R/hazard_api.R` ~L404–415), so
they surface at fit time — this is the "loud" goal. The pre-existing
`"Hessian not invertible; standard errors unavailable"` string is reused on
the non-invertible fallback path and stays muffled by the existing handler.
`summary()` additionally shows a note for users who fit elsewhere.

---

## Task 1: `.hzr_safe_solve()` helper

**Files:**
- Create: `R/hessian-invert.R`
- Test: `tests/testthat/test-hessian-invert.R`

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-hessian-invert.R`:

```r
test_that("well-conditioned PD Hessian inverts cleanly with diagnostics", {
  H <- matrix(c(4, 1, 1, 3), 2, 2)
  res <- expect_silent(.hzr_safe_solve(H))
  expect_equal(res$vcov, solve(H), tolerance = 1e-10)
  expect_true(res$pd)
  expect_true(is.finite(res$rcond) && res$rcond > 1e-3)
})

test_that("asymmetric input is symmetrized before inversion", {
  H_sym  <- matrix(c(4, 1, 1, 3), 2, 2)
  H_asym <- matrix(c(4, 0, 2, 3), 2, 2)  # (H+t(H))/2 == H_sym
  expect_equal(.hzr_safe_solve(H_asym)$vcov, solve(H_sym), tolerance = 1e-10)
})

test_that("ill-conditioned Hessian warns by name but still returns numbers", {
  H <- matrix(c(1, 1, 1, 1 + 1e-12), 2, 2)
  expect_warning(res <- .hzr_safe_solve(H), "ill-conditioned")
  expect_true(is.matrix(res$vcov))
  expect_true(is.finite(res$rcond))
})

test_that("non-PD Hessian falls back to solve and reports pd = FALSE", {
  # Indefinite (eigenvalues +/-): chol() fails, solve() succeeds.
  H <- matrix(c(1, 2, 2, 1), 2, 2)
  res <- suppressWarnings(.hzr_safe_solve(H))
  expect_false(res$pd)
  expect_true(is.matrix(res$vcov))
})

test_that("non-positive variance diagonals are NA'd with a warning", {
  # Indefinite matrix whose inverse has a negative diagonal entry.
  H <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_warning(res <- .hzr_safe_solve(H), "[Nn]on-positive variance")
  expect_true(any(is.na(diag(res$vcov))))
})

test_that("non-finite Hessian returns NA vcov with a warning", {
  H <- matrix(c(1, NA, NA, 1), 2, 2)
  expect_warning(res <- .hzr_safe_solve(H), "non-finite")
  expect_true(all(is.na(res$vcov)))
  expect_true(is.na(res$pd))
})

test_that("scalar / non-matrix input is handled as non-finite", {
  expect_warning(res <- .hzr_safe_solve(NA), "non-finite")
  expect_true(all(is.na(res$vcov)))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hessian-invert.R")'`
Expected: FAIL — `could not find function ".hzr_safe_solve"`.

- [ ] **Step 3: Write the implementation**

Create `R/hessian-invert.R`:

```r
#' @keywords internal
NULL

# hessian-invert.R -- Numerically stable Hessian inversion with diagnostics.

# Reciprocal-condition warning threshold (~1.5e-8). Shared by .hzr_safe_solve()
# and the summary() diagnostic note so they never drift apart.
.hzr_rcond_tol <- .Machine$double.eps^0.5

#' Stable Hessian inversion with conditioning diagnostics
#'
#' Inverts a negative-log-likelihood Hessian into a variance-covariance
#' matrix, hardening against the ill-conditioning that arises at high
#' parameter counts.  Symmetrizes the input, checks the reciprocal condition
#' number, inverts via Cholesky (with a \code{solve()} fallback for non-PD
#' Hessians), and guards non-positive variances.  Emits a named warning for
#' every degenerate path.
#'
#' @param H Square numeric Hessian of the negative log-likelihood.
#' @param tol Reciprocal-condition warning threshold.
#' @return A list with:
#'   \code{vcov} (the variance-covariance matrix, or \code{NA} on failure;
#'   diagonals with non-positive variance, and their rows/cols, set to
#'   \code{NA}); \code{rcond} (reciprocal condition number of the symmetrized
#'   Hessian, \code{NA} if unavailable); \code{pd} (\code{TRUE} if the Hessian
#'   was positive-definite, \code{FALSE} if inverted via fallback, \code{NA}
#'   if not invertible).
#' @noRd
.hzr_safe_solve <- function(H, tol = .hzr_rcond_tol) {
  # (1) Non-finite / non-matrix guard
  if (is.null(H) || !is.matrix(H) || anyNA(H) || any(!is.finite(H))) {
    warning("Hessian contains non-finite entries; standard errors unavailable")
    return(list(vcov = NA, rcond = NA_real_, pd = NA))
  }

  # (2) Symmetrize (numDeriv Hessians are only symmetric to Richardson tol)
  H <- (H + t(H)) / 2

  # (3) Conditioning check
  rc <- tryCatch(rcond(H), error = function(e) NA_real_)
  if (is.na(rc) || rc < tol) {
    warning(sprintf(
      "Hessian is ill-conditioned (rcond = %.3g); standard errors may be unreliable",
      rc
    ))
  }

  # (4) Stable inversion: Cholesky (PD) with solve() fallback (non-PD)
  ch <- tryCatch(chol(H), error = function(e) NULL)
  if (!is.null(ch)) {
    pd <- TRUE
    vcov <- chol2inv(ch)
  } else {
    pd <- FALSE
    vcov <- tryCatch(solve(H), error = function(e) NULL)
    if (is.null(vcov)) {
      warning("Hessian not invertible; standard errors unavailable")
      return(list(vcov = NA, rcond = rc, pd = NA))
    }
  }
  dimnames(vcov) <- dimnames(H)  # chol2inv() drops names

  # (5) Non-positive-variance guard
  d <- diag(vcov)
  bad <- !is.finite(d) | d < 0
  if (any(bad)) {
    warning("Non-positive variance estimates; the optimum may not be a proper maximum")
    vcov[bad, ] <- NA_real_
    vcov[, bad] <- NA_real_
  }

  list(vcov = vcov, rcond = rc, pd = pd)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hessian-invert.R")'`
Expected: PASS — all expectations pass, 0 failures.

- [ ] **Step 5: Commit**

```bash
git add R/hessian-invert.R tests/testthat/test-hessian-invert.R
git commit -m "feat: add .hzr_safe_solve() hardened Hessian inversion (Layer 1)"
```

---

## Task 2: Wire `.hzr_safe_solve()` into the chokepoint

**Files:**
- Modify: `R/optimizer.R` (~L129–162)
- Test: `tests/testthat/test-hessian-invert.R` (append)

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-hessian-invert.R`:

```r
test_that(".hzr_optim_generic returns rcond and pd diagnostics", {
  set.seed(1)
  n <- 200L
  time <- rexp(n, rate = 0.5)
  status <- rep(1L, n)
  res <- .hzr_optim_generic(
    logl_fn = .hzr_logl_exponential,
    gradient_fn = .hzr_gradient_exponential,
    time = time, status = status,
    x = NULL, theta_start = c(log_rate = 0),
    control = list(maxit = 200, reltol = 1e-8, abstol = 1e-8),
    use_bounds = FALSE
  )
  expect_true(is.finite(res$rcond))
  expect_true(isTRUE(res$pd))
  expect_true(is.matrix(res$vcov) && all(is.finite(diag(res$vcov))))
})
```

(Note: `theta_start` name/length follows the exponential parameterization;
adjust the name only if `.hzr_logl_exponential` expects a different label —
the value `0` and length 1 are correct for a covariate-free exponential fit.)

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hessian-invert.R")'`
Expected: FAIL — `res$rcond` is `NULL` (not yet returned by the chokepoint).

- [ ] **Step 3: Replace the inversion block in `.hzr_optim_generic()`**

In `R/optimizer.R`, replace the existing post-fit Hessian + vcov block
(currently lines ~129–162, from the `# Post-fit Hessian for standard errors`
comment through the final `list(...)` return) with:

```r
  # Post-fit Hessian for standard errors
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

  # Hardened inversion + conditioning diagnostics (Layer 1).
  inv <- if (is.matrix(hess_result)) {
    .hzr_safe_solve(hess_result)
  } else {
    list(vcov = NA, rcond = NA_real_, pd = NA)
  }

  list(
    par = result$par,
    value = -result$value,
    convergence = result$convergence,
    counts = result$counts,
    message = result$message,
    hessian = hess_result,
    vcov = inv$vcov,
    rcond = inv$rcond,
    pd = inv$pd
  )
}
```

(The guard `if (is.matrix(hess_result))` avoids a spurious "non-finite"
warning when `numDeriv` is unavailable and `hess_result` is scalar `NA`.)

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hessian-invert.R")'`
Expected: PASS — all expectations including the new diagnostics test.

- [ ] **Step 5: Run the full optimizer-touching suites for regression**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat", filter = "weibull|exponential|loglogistic|lognormal|multiphase|wald")'`
Expected: 0 failures (vcov values unchanged for well-conditioned fits —
`chol2inv` and `solve` agree to numerical tolerance).

- [ ] **Step 6: Commit**

```bash
git add R/optimizer.R tests/testthat/test-hessian-invert.R
git commit -m "feat: route .hzr_optim_generic vcov through .hzr_safe_solve; emit rcond/pd"
```

---

## Task 3: Thread `rcond`/`pd` onto the fit object

**Files:**
- Modify: `R/hazard_api.R` (~L444 multiphase branch, ~L474 single-dist branch)
- Test: `tests/testthat/test-hessian-stability.R` (create)

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-hessian-stability.R`:

```r
test_that("fit object carries rcond and pd diagnostics (weibull)", {
  set.seed(2)
  n <- 250L
  dat <- data.frame(t = rweibull(n, shape = 1.4, scale = 3),
                    d = rep(1L, n))
  fit <- hazard(survival::Surv(t, d) ~ 1, data = dat, dist = "weibull",
                fit = TRUE)
  expect_true(is.finite(fit$fit$rcond))
  expect_true(isTRUE(fit$fit$pd))
})

test_that("fit object carries rcond and pd diagnostics (multiphase)", {
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                    fixed = "m"),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 3, maxit = 500, conserve = TRUE)
  )
  expect_true(is.finite(fit$fit$rcond))
  expect_false(is.null(fit$fit$pd))
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hessian-stability.R")'`
Expected: FAIL — `fit$fit$rcond` is `NULL`.

- [ ] **Step 3: Copy diagnostics into `fit_state` (both branches)**

In `R/hazard_api.R`, in the **multiphase** branch, immediately after the line
`fit_state$vcov <- optim_result$vcov` (~L444), add:

```r
    fit_state$rcond <- optim_result$rcond
    fit_state$pd <- optim_result$pd
```

In the **single-distribution** branch, immediately after the line
`fit_state$vcov <- optim_result$vcov` (~L474), add the identical two lines:

```r
    fit_state$rcond <- optim_result$rcond
    fit_state$pd <- optim_result$pd
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hessian-stability.R")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/hazard_api.R tests/testthat/test-hessian-stability.R
git commit -m "feat: thread rcond/pd diagnostics onto the fit object"
```

---

## Task 4: Flagged diagnostic note in `summary()`

**Files:**
- Modify: `R/hazard_api.R` (`summary.hazard` out-list ~L1121; `print.summary.hazard` ~L1180)
- Test: `tests/testthat/test-hessian-stability.R` (append)

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-hessian-stability.R`:

```r
test_that("summary prints no Hessian note for a clean fit", {
  set.seed(3)
  dat <- data.frame(t = rweibull(250, 1.4, 3), d = rep(1L, 250))
  fit <- hazard(survival::Surv(t, d) ~ 1, data = dat, dist = "weibull",
                fit = TRUE)
  out <- capture.output(print(summary(fit)))
  expect_false(any(grepl("ill-conditioned|not positive-definite", out)))
})

test_that("summary prints a note when the fit is flagged ill-conditioned", {
  s <- summary(structure(
    list(spec = list(dist = "weibull", phases = NULL),
         engine = "test",
         fit = list(theta = c(a = 1), vcov = matrix(1, 1, 1),
                    converged = TRUE, objective = -1,
                    rcond = 1e-12, pd = TRUE),
         data = list(time = 1:5, x = NULL),
         call = quote(hazard())),
    class = "hazard"))
  out <- capture.output(print(s))
  expect_true(any(grepl("ill-conditioned", out)))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hessian-stability.R")'`
Expected: FAIL — the flagged-fit test finds no "ill-conditioned" line.

- [ ] **Step 3: Store diagnostics in the summary object**

In `R/hazard_api.R`, in `summary.hazard()`, add two fields to the `out <-
list(...)` (after `has_vcov = ...`, ~L1132):

```r
    rcond = object$fit$rcond,
    pd = object$fit$pd,
```

- [ ] **Step 4: Print the flagged note**

In `print.summary.hazard()`, immediately after the log-lik block (the
`if (!is.null(x$log_lik) ...)` block ending ~L1182), add:

```r
  if (!is.null(x$rcond) && !is.na(x$rcond) && x$rcond < .hzr_rcond_tol) {
    cat("  Note: Hessian ill-conditioned (rcond = ",
        format(x$rcond, digits = 3),
        "); standard errors may be unreliable.\n", sep = "")
  }
  if (!is.null(x$pd) && !is.na(x$pd) && !isTRUE(x$pd)) {
    cat("  Note: Hessian not positive-definite at the optimum; ",
        "standard errors may be unreliable.\n", sep = "")
  }
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hessian-stability.R")'`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add R/hazard_api.R tests/testthat/test-hessian-stability.R
git commit -m "feat: summary() prints a note when the Hessian is flagged"
```

---

## Task 5: 13-parameter multiphase anchor regression

**Files:**
- Test: `tests/testthat/test-hessian-stability.R` (append)

This is the row's named borderline case. It pins the post-Layer-1 behavior:
the high-dimensional fit converges and produces finite, positive SEs with a
recoverable conditioning diagnostic.

- [ ] **Step 1: Write the anchor test**

Append to `tests/testthat/test-hessian-stability.R`:

```r
test_that("13-parameter multiphase deciles fit is stable (anchor)", {
  skip_on_cran()
  data(avc, package = "TemporalHazard")
  avc <- na.omit(avc)

  # A high-dimensional multiphase fit with covariates in both phases,
  # exercising the 12+-parameter inversion path named in DEVELOPMENT-PLAN
  # §7c ("hm.death.AVC.deciles").
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", formula = ~ age + female + inc_surg,
                           t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
      constant = hzr_phase("constant", formula = ~ age + female + inc_surg)
    ),
    fit = TRUE, control = list(n_starts = 3, maxit = 800, conserve = TRUE)
  )

  expect_true(fit$fit$converged)

  free <- which(is.finite(diag(fit$fit$vcov)))
  expect_gte(length(free), 12L)               # genuinely high-dimensional
  se <- sqrt(diag(fit$fit$vcov)[free])
  expect_true(all(is.finite(se)))             # no NaN/NA SEs on free params
  expect_true(all(se > 0))                    # positive variances
  expect_true(is.finite(fit$fit$rcond))       # conditioning was measured
})
```

(Note: the exact covariate names `age`, `female`, `inc_surg` are columns of
the shipped `avc` dataset used elsewhere in the suite, e.g.
`test-phase-specific-covariates.R`. Verify column names with
`names(avc)` during Step 2; substitute the dataset's actual covariate
columns if any differ, keeping the free-parameter count ≥ 12.)

- [ ] **Step 2: Run the anchor test**

Run: `Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-hessian-stability.R")'`
Expected: PASS. If it FAILS on `expect_gte(length(free), 12L)`, the chosen
covariates yielded < 12 free params — add covariate columns until ≥ 12. If
it FAILS on `pd`/`rcond`, that is a genuine finding: the fit is legitimately
ill-conditioned (not a numDeriv artifact); per the design's open question,
soften the anchor to `expect_warning(..., "ill-conditioned")` / assert the
summary note instead of asserting cleanliness, and record the finding in the
PR description.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-hessian-stability.R
git commit -m "test: 13-parameter multiphase Hessian-stability anchor"
```

---

## Task 6: NEWS entry + full check

**Files:**
- Modify: `NEWS.md`

- [ ] **Step 1: Add the NEWS entry**

In `NEWS.md`, under the existing `# TemporalHazard 1.1.0.9000 (development
version)` heading, in the `## Bug fixes` (or a new `## Improvements`)
subsection, add:

```markdown
* **Hardened Hessian inversion for standard errors (Phase 7c).**
  Post-fit variance-covariance estimation now symmetrizes the Hessian,
  checks its reciprocal condition number, inverts via Cholesky with a
  `solve()` fallback for non-positive-definite Hessians, and guards
  non-positive variances instead of silently emitting `NaN` standard
  errors. Ill-conditioned, non-positive-definite, and non-finite Hessians
  now raise specific, named warnings, and fits carry `rcond` / `pd`
  diagnostics that `summary()` surfaces as a note when a fit is flagged.
  This closes the "12+-parameter Hessian stability" hardening item for the
  inversion layer; analytic Hessians (more accurate standard errors) follow
  in subsequent releases.
```

- [ ] **Step 2: Run the full test suite**

Run: `Rscript -e 'devtools::test()'`
Expected: 0 failures across the whole suite.

- [ ] **Step 3: Run R CMD check (lint/Rd/usage sanity)**

Run: `Rscript -e 'devtools::check(args = c("--no-manual"), quiet = TRUE)'`
Expected: 0 errors / 0 warnings; the pre-existing single NOTE only. (Confirm
no new `object_usage_linter` notes from the new internal `.hzr_safe_solve` /
`.hzr_rcond_tol` symbols.)

- [ ] **Step 4: Commit**

```bash
git add NEWS.md
git commit -m "docs: NEWS entry for hardened Hessian inversion (Layer 1)"
```

---

## Self-Review Notes

- **Spec coverage (Layer 1 portion of HESSIAN-STABILITY-DESIGN.md):**
  symmetrize (T1 §2), finiteness guard (T1 §1), rcond check + warn-but-return
  policy (T1 §3), Cholesky-with-fallback (T1 §4), negative-variance guard
  (T1 §5), `rcond`/`pd` on fit (T3), summary note (T4), Layer-1 unit tests
  (T1), 13-param anchor (T5). Analytic-Hessian items (cross-check / scale-
  invariance / `hessian_fn` plumbing) are deferred to PRs 2–6 by design and
  intentionally absent here.
- **Threshold DRY:** `.hzr_rcond_tol` is defined once in `R/hessian-invert.R`
  and reused by both `.hzr_safe_solve()` and `print.summary.hazard()`.
- **No production regression:** for well-conditioned fits, `chol2inv(chol(H))`
  equals `solve(H)` to numerical tolerance, so existing vcov/SE-asserting
  tests (e.g. `test-wald.R`, `test-sas-parity.R`) should not move; Task 2
  Step 5 verifies this before the change propagates further.
