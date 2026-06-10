# Weighted Multiphase + Covariates — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close roadmap item 7c ("weighted multiphase + covariates") by adding integer-weight duplication-parity tests that exercise a covariate entering a phase at the likelihood, analytic-gradient, and end-to-end fit levels (CoE off and on).

**Architecture:** Test-only. No `R/` changes unless a test surfaces a bug. Reuse the existing duplication invariant `L(theta; D, w) == L(theta; expand(D, w), 1)`, extending it from intercept-only to covariate-bearing phases via `hzr_phase(formula = ~ x)` and explicit `covariate_counts`/`x_list` for the direct-LL calls.

**Tech Stack:** R, testthat, numDeriv. Internal functions `TemporalHazard:::.hzr_logl_multiphase()`, `:::.hzr_gradient_multiphase()`; public `hazard(dist = "multiphase", ...)`.

**Critical execution rule:** These tests assert *existing* code is correct, so each should PASS on first run. **A FAIL is a bug discovery, not a test to massage** — stop and invoke superpowers:systematic-debugging. Never relax the invariant or loosen a tolerance to make a red test green.

**Theta layout reference** (early `cdf` with `fixed = "shapes"` + one covariate, then `constant` + one covariate):
`[log_mu_early, log_t_half, nu, m, beta_early, log_mu_constant, beta_constant]` — per phase the covariate `beta`s come *after* that phase's shape block (see `.hzr_unpack_phase_theta`, `R/likelihood-multiphase.R:122`).

---

### Task 1: Likelihood duplication parity with a covariate

**Files:**
- Modify: `tests/testthat/test-weights.R` (append after the existing multiphase-likelihood test, ~line 203)

- [ ] **Step 1: Add a shared fixture helper + the likelihood test**

Insert after the `"integer weights match duplication in multiphase likelihood"` block (line 203):

```r
# ---------------------------------------------------------------------------
# Multiphase WITH covariates: weighted == duplicated (roadmap 7c)
# ---------------------------------------------------------------------------

# Shared fixture: multiphase data with a covariate in BOTH phases, integer
# weights, and the row-duplication index. covariate enters early + constant.
make_mp_cov <- function(seed = 41, n = 40) {
  set.seed(seed)
  t      <- runif(n, 0.05, 3)
  status <- rbinom(n, 1, 0.4)
  x      <- rnorm(n)
  w      <- sample(1:3, n, replace = TRUE)
  idx    <- rep(seq_len(n), times = w)
  list(t = t, status = status, x = x, w = w, idx = idx, n = n)
}

# Phases + theta shared by the direct-LL and gradient tests.
# Layout: [log_mu_e, log_thalf, nu, m, beta_e, log_mu_c, beta_c]
mp_cov_phases <- function() {
  list(
    early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
}
mp_cov_theta   <- c(-4, log(0.3), 1, 1, 0.5, -3, -0.3)
mp_cov_counts  <- c(early = 1L, constant = 1L)

test_that("integer weights match duplication in multiphase LL WITH covariates", {
  skip_on_cran()
  d <- make_mp_cov()
  x_list   <- list(early    = matrix(d$x, ncol = 1),
                   constant = matrix(d$x, ncol = 1))
  x_list_d <- list(early    = matrix(d$x[d$idx], ncol = 1),
                   constant = matrix(d$x[d$idx], ncol = 1))

  ll_w <- TemporalHazard:::.hzr_logl_multiphase(
    theta = mp_cov_theta, time = d$t, status = d$status,
    phases = mp_cov_phases(), covariate_counts = mp_cov_counts,
    x_list = x_list, weights = d$w
  )
  ll_dup <- TemporalHazard:::.hzr_logl_multiphase(
    theta = mp_cov_theta, time = d$t[d$idx], status = d$status[d$idx],
    phases = mp_cov_phases(), covariate_counts = mp_cov_counts,
    x_list = x_list_d, weights = rep(1, length(d$idx))
  )
  expect_equal(ll_w, ll_dup, tolerance = 1e-10)
})
```

- [ ] **Step 2: Run the test**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-weights.R", filter = NULL)'`
or focused: `Rscript -e 'devtools::test(filter = "weights")'`
Expected: PASS. If FAIL → STOP, invoke systematic-debugging (the weighted covariate LL is wrong).

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-weights.R
git commit -m "test(weights): multiphase LL duplication parity with covariates (7c)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: Analytic gradient vs numerical, with a covariate

**Files:**
- Modify: `tests/testthat/test-weights.R` (append after Task 1's test)

- [ ] **Step 1: Add the gradient test**

```r
test_that("weighted multiphase analytic gradient matches numerical WITH covariates", {
  skip_if_not_installed("numDeriv")
  d <- make_mp_cov(seed = 42)
  x_list <- list(early    = matrix(d$x, ncol = 1),
                 constant = matrix(d$x, ncol = 1))

  num_g <- numDeriv::grad(
    function(th) {
      TemporalHazard:::.hzr_logl_multiphase(
        theta = th, time = d$t, status = d$status,
        phases = mp_cov_phases(), covariate_counts = mp_cov_counts,
        x_list = x_list, weights = d$w
      )
    },
    mp_cov_theta
  )
  ana_g <- TemporalHazard:::.hzr_gradient_multiphase(
    theta = mp_cov_theta, time = d$t, status = d$status, weights = d$w,
    phases = mp_cov_phases(), covariate_counts = mp_cov_counts, x_list = x_list
  )
  expect_equal(as.numeric(ana_g), num_g, tolerance = 1e-4)
})
```

- [ ] **Step 2: Run the test**

Run: `Rscript -e 'devtools::test(filter = "weights")'`
Expected: PASS. A FAIL means the analytic beta-score under weights is wrong → systematic-debugging. (The beta-score is the column never exercised by the intercept-only gradient test.)

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-weights.R
git commit -m "test(weights): weighted multiphase gradient vs numerical with covariates (7c)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: End-to-end weighted fit == duplicated fit (CoE OFF)

**Files:**
- Modify: `tests/testthat/test-weights.R` (append after Task 2's test)

- [ ] **Step 1: Add the end-to-end fit test (conserve = FALSE)**

Covariate enters both phases via `hzr_phase(formula = ~ x)`. The duplicated
data frame is built with the existing `expand_rows()` helper.

```r
test_that("weighted multiphase fit matches duplicated-row fit, covariates, CoE off", {
  skip_on_cran()
  d  <- make_mp_cov(seed = 43)
  df <- data.frame(time = d$t, status = d$status, x = d$x)
  df_exp <- expand_rows(df, d$w)

  mk_fit <- function(data, weights) {
    hazard(
      survival::Surv(time, status) ~ 1,
      data = data, dist = "multiphase",
      phases = list(
        early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                             fixed = "shapes", formula = ~ x),
        constant = hzr_phase("constant", formula = ~ x)
      ),
      weights = weights, fit = TRUE,
      control = list(conserve = FALSE, n_starts = 1L, maxit = 500L)
    )
  }

  fit_w   <- mk_fit(df,     d$w)
  fit_dup <- mk_fit(df_exp, NULL)

  expect_equal(coef(fit_w), coef(fit_dup), tolerance = 1e-3)
  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-4)
})
```

- [ ] **Step 2: Run the test**

Run: `Rscript -e 'devtools::test(filter = "weights")'`
Expected: PASS. A FAIL here (when Tasks 1–2 passed) points at the optimizer/Hessian path under weights+covariates → systematic-debugging.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-weights.R
git commit -m "test(weights): end-to-end weighted multiphase+covariate fit parity, CoE off (7c)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4: End-to-end weighted fit == duplicated fit (CoE ON)

**Files:**
- Modify: `tests/testthat/test-conservation-of-events.R` (append at end of file)

- [ ] **Step 1: Add the CoE-on covariate fit-parity test**

This is the new combination the roadmap flags: weighted Conservation-of-Events
*with* a covariate. Identical to Task 3 but `conserve = TRUE`.

```r
test_that("weighted multiphase fit matches duplicated-row fit, covariates, CoE on", {
  skip_on_cran()
  set.seed(43)
  n      <- 40
  t      <- runif(n, 0.05, 3)
  status <- rbinom(n, 1, 0.4)
  x      <- rnorm(n)
  w      <- sample(1:3, n, replace = TRUE)
  df     <- data.frame(time = t, status = status, x = x)
  df_exp <- df[rep(seq_len(n), times = w), , drop = FALSE]

  mk_fit <- function(data, weights) {
    hazard(
      survival::Surv(time, status) ~ 1,
      data = data, dist = "multiphase",
      phases = list(
        early    = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                             fixed = "shapes", formula = ~ x),
        constant = hzr_phase("constant", formula = ~ x)
      ),
      weights = weights, fit = TRUE,
      control = list(conserve = TRUE, n_starts = 1L, maxit = 500L)
    )
  }

  fit_w   <- mk_fit(df,     w)
  fit_dup <- mk_fit(df_exp, NULL)

  expect_equal(coef(fit_w), coef(fit_dup), tolerance = 1e-3)
  expect_equal(fit_w$fit$objective, fit_dup$fit$objective, tolerance = 1e-4)
})
```

- [ ] **Step 2: Run the test**

Run: `Rscript -e 'devtools::test(filter = "conservation-of-events")'`
Expected: PASS. A FAIL isolates the weighted-CoE-with-covariate interaction (`.hzr_conserve_events()` / `.hzr_select_fixmu_phase()` scaling under a covariate-bearing phase) → systematic-debugging.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-conservation-of-events.R
git commit -m "test(coe): weighted multiphase+covariate fit parity, CoE on (7c)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 5: Full-suite verification + roadmap close-out

**Files:**
- Modify: `inst/dev/DEVELOPMENT-PLAN.md` (the 7c "Weighted multiphase + covariates" row, line ~530)

- [ ] **Step 1: Run the full test suite**

Run: `Rscript -e 'devtools::test()'`
Expected: 0 failures. Confirm the four new tests are counted and green, and nothing else regressed.

- [ ] **Step 2: Flip the 7c roadmap row to covered**

In the 7c table, change the "Weighted multiphase + covariates" row Status from
`No test` to `✅ Tested (LL + gradient + fit, CoE on/off)`.

- [ ] **Step 3: Commit**

```bash
git add inst/dev/DEVELOPMENT-PLAN.md
git commit -m "docs(roadmap): mark 7c weighted multiphase+covariates as tested

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

- [ ] **Step 4: Open the PR**

```bash
git push -u origin test/weighted-multiphase-covariates
gh pr create --base main --title "test: weighted multiphase + covariates coverage (7c)" \
  --body "Closes roadmap 7c gap: integer-weight duplication parity for multiphase fits with covariates in a phase — likelihood, analytic gradient, and end-to-end fit (CoE off and on). Test-only; no R/ changes."
```

---

## Notes for the implementer

- If `devtools` is unavailable, substitute `testthat::test_local()` or
  `R CMD INSTALL . && Rscript -e 'testthat::test_dir("tests/testthat")'`.
- The `mp_cov_*` helpers (Task 1) are defined once and reused by Tasks 1–3;
  do not redefine them. Task 4 lives in a different file and is intentionally
  self-contained (no cross-file helper sharing in testthat).
- Tolerances (`1e-10` LL, `1e-4` gradient/objective, `1e-3` coef) mirror the
  surrounding intercept-only tests — do not tighten or loosen them.
