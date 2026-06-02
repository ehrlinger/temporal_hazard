# Per-phase CLs for predict(decompose=TRUE, se.fit=TRUE) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Unblock `predict(type="cumulative_hazard", decompose=TRUE, se.fit=TRUE)` for multiphase models, returning long-format per-phase + total cumulative-hazard delta-method confidence limits.

**Architecture:** Reuse the existing per-phase Jacobian blocks already built by `.hzr_predict_jacobian_multiphase()` and the existing delta-method sandwich. Add a `per_phase` mode to that Jacobian builder, extract the shared vcov free-parameter restriction into a helper, add a decomposed se.fit driver, and dispatch to it from `predict.hazard()` instead of the current hard error.

**Tech Stack:** R, testthat, numDeriv, roxygen2. Internal functions in `R/predict-cl.R`; public `predict.hazard()` in `R/hazard_api.R`.

**Spec:** `inst/dev/PREDICT-DECOMPOSE-SE-DESIGN.md`

**Key math:** `H(t|x) = sum_j mu_j(x) Phi_j(t)`. Each phase's contribution `H_j = mu_j Phi_j` depends only on phase j's parameters, so its Jacobian block has nonzero columns only for phase j. Per-phase SE `= sqrt(diag(J_j V J_j^T))`; the total uses the summed Jacobian. Per-phase CLs do NOT sum to the total CL (cross-phase covariance lives only in the total) — this is correct and documented.

**File map:**
- `R/predict-cl.R` — add `per_phase` arg (Task 1), extract `.hzr_free_vcov()` (Task 2), add `.hzr_predict_with_se_decomposed()` (Task 3).
- `R/hazard_api.R` — dispatch + roxygen (Task 4, Task 5).
- `tests/testthat/test-predict-cl.R` — tests across tasks.
- `NEWS.md`, `man/*.Rd` — Task 5.

---

### Task 1: `per_phase` mode on the multiphase Jacobian builder

**Files:**
- Modify: `R/predict-cl.R` (`.hzr_predict_jacobian_multiphase`, ~lines 143-208)
- Test: `tests/testthat/test-predict-cl.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-predict-cl.R`:

```r
# ---------------------------------------------------------------------------
# per_phase Jacobian mode (decompose + se.fit support)
# ---------------------------------------------------------------------------

test_that("multiphase jacobian per_phase=TRUE returns blocks that sum to the total", {
  skip_on_cran()
  t_new <- c(0.3, 0.8, 1.5, 2.2)
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  cov_counts <- c(early = 0L, constant = 0L)
  x_list <- list(early = NULL, constant = NULL)
  theta <- c(-3, log(0.3), 1, 1, -2)
  p <- length(theta)

  J_total <- TemporalHazard:::.hzr_predict_jacobian_multiphase(
    theta, t_new, phases, cov_counts, x_list, p
  )
  J_list <- TemporalHazard:::.hzr_predict_jacobian_multiphase(
    theta, t_new, phases, cov_counts, x_list, p, per_phase = TRUE
  )

  expect_type(J_list, "list")
  expect_named(J_list, c("early", "constant"))
  expect_equal(dim(J_list$early), c(length(t_new), p))
  # Blocks sum to the aggregate Jacobian.
  expect_equal(Reduce(`+`, J_list), J_total, tolerance = 1e-12)
  # Each block touches only its own phase's columns:
  # early occupies cols 1:4 (log_mu, log_thalf, nu, m), constant col 5 (log_mu).
  expect_true(all(J_list$early[, 5] == 0))
  expect_true(all(J_list$constant[, 1:4] == 0))
})
```

- [ ] **Step 2: Run it; expect FAIL**

Run: `Rscript -e 'suppressMessages(devtools::load_all(quiet=TRUE)); testthat::test_file("tests/testthat/test-predict-cl.R")'`
Expected: FAIL — `per_phase` is an unused argument / list not returned.

- [ ] **Step 3: Implement `per_phase`**

Replace the body of `.hzr_predict_jacobian_multiphase` (signature + internals) in `R/predict-cl.R` with:

```r
.hzr_predict_jacobian_multiphase <- function(theta, time, phases,
                                              covariate_counts, x_list, p,
                                              per_phase = FALSE) {
  n <- length(time)
  theta_split <- .hzr_split_theta(theta, phases, covariate_counts)
  # One n x p matrix per phase; each phase fills only its own columns.
  J_list <- stats::setNames(
    lapply(seq_along(phases), function(i) matrix(0, nrow = n, ncol = p)),
    names(phases)
  )

  pos <- 1L
  for (nm in names(phases)) {
    ph <- phases[[nm]]
    pars <- .hzr_unpack_phase_theta(theta_split[[nm]], ph)
    Jp <- J_list[[nm]]

    # mu_j(x_i)
    if (length(pars$beta) > 0L && !is.null(x_list[[nm]])) {
      eta_j <- pars$log_mu + as.numeric(x_list[[nm]] %*% pars$beta)
    } else {
      eta_j <- rep(pars$log_mu, n)
    }
    mu_j <- exp(eta_j)

    # Phi_j and shape derivatives
    if (ph$type == "constant") {
      Phi_j <- time
      Jp[, pos] <- mu_j * Phi_j
      pos <- pos + 1L
    } else if (ph$type == "g3") {
      tau_j <- exp(pars$log_tau)
      pd <- .hzr_g3_phase_derivatives(time, tau = tau_j,
                                        gamma = pars$gamma,
                                        alpha = pars$alpha,
                                        eta = pars$eta)
      Phi_j <- pd$Phi
      Jp[, pos] <- mu_j * Phi_j
      Jp[, pos + 1L] <- mu_j * pd$dPhi_dlog_tau
      Jp[, pos + 2L] <- mu_j * pd$dPhi_dgamma
      Jp[, pos + 3L] <- mu_j * pd$dPhi_dalpha
      Jp[, pos + 4L] <- mu_j * pd$dPhi_deta
      pos <- pos + 5L
    } else {
      t_half_j <- exp(pars$log_t_half)
      pd <- .hzr_phase_derivatives(time, t_half = t_half_j,
                                    nu = pars$nu, m = pars$m,
                                    type = ph$type)
      Phi_j <- pd$Phi
      Jp[, pos] <- mu_j * Phi_j
      Jp[, pos + 1L] <- mu_j * pd$dPhi_dlog_thalf
      Jp[, pos + 2L] <- mu_j * pd$dPhi_dnu
      Jp[, pos + 3L] <- mu_j * pd$dPhi_dm
      pos <- pos + 4L
    }

    # Covariate betas
    n_beta <- covariate_counts[[nm]]
    if (n_beta > 0L && !is.null(x_list[[nm]])) {
      x_phase <- x_list[[nm]]
      for (k in seq_len(n_beta)) {
        Jp[, pos] <- x_phase[, k] * mu_j * Phi_j
        pos <- pos + 1L
      }
    } else if (n_beta > 0L) {
      pos <- pos + n_beta
    }

    J_list[[nm]] <- Jp
  }

  if (per_phase) {
    return(J_list)
  }
  Reduce(`+`, J_list)
}
```

Also update the roxygen `@return` line above the function to:

```r
#' @param per_phase Logical; if `TRUE`, return a named list of per-phase
#'   `n x p` Jacobians (each with only that phase's columns nonzero) instead
#'   of their sum.
#' @return Numeric `n x p` Jacobian of `H(t|x)` (default), or a named list of
#'   per-phase `n x p` Jacobians when `per_phase = TRUE`.
```

- [ ] **Step 4: Run it; expect PASS**

Run: `Rscript -e 'suppressMessages(devtools::load_all(quiet=TRUE)); testthat::test_file("tests/testthat/test-predict-cl.R")'`
Expected: PASS, including the pre-existing "Multiphase analytic jacobian matches numDeriv" test (per_phase=FALSE path unchanged).

- [ ] **Step 5: Commit**

```bash
git add R/predict-cl.R tests/testthat/test-predict-cl.R
git commit -m "feat(predict-cl): per_phase mode on multiphase Jacobian builder

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: Extract shared vcov free-parameter helper

**Files:**
- Modify: `R/predict-cl.R` (`.hzr_predict_with_se`, vcov block ~lines 356-401; add `.hzr_free_vcov()`)
- Test: `tests/testthat/test-predict-cl.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-predict-cl.R`:

```r
test_that(".hzr_free_vcov restricts to finite-diagonal free parameters", {
  # All-free 2x2.
  V <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
  res <- TemporalHazard:::.hzr_free_vcov(V, p = 2L)
  expect_equal(res$free_idx, 1:2)
  expect_equal(res$vcov_use, V)

  # One fixed parameter (NA row/col) -> restricted to the free submatrix.
  Vf <- matrix(NA_real_, 3, 3)
  Vf[c(1, 3), c(1, 3)] <- c(0.04, 0.00, 0.00, 0.09)
  res2 <- TemporalHazard:::.hzr_free_vcov(Vf, p = 3L)
  expect_equal(res2$free_idx, c(1L, 3L))
  expect_equal(dim(res2$vcov_use), c(2L, 2L))

  # Unusable vcov -> NULL with a warning.
  expect_warning(
    bad <- TemporalHazard:::.hzr_free_vcov(NULL, p = 3L),
    "unavailable"
  )
  expect_null(bad)
})
```

- [ ] **Step 2: Run it; expect FAIL**

Run: `Rscript -e 'suppressMessages(devtools::load_all(quiet=TRUE)); testthat::test_file("tests/testthat/test-predict-cl.R")'`
Expected: FAIL — `.hzr_free_vcov` not found.

- [ ] **Step 3: Add the helper and refactor `.hzr_predict_with_se` to use it**

Add this function in `R/predict-cl.R` directly above `.hzr_predict_with_se` (i.e. before the `#' Compute point predictions + delta-method CLs` block):

```r
#' Free-parameter vcov submatrix for the delta-method sandwich
#'
#' Fixed parameters (from `fixed = "shapes"` or CoE) leave NA rows/cols in the
#' expanded vcov. Treat them as known-with-zero-variance: restrict the sandwich
#' to the free submatrix. Returns `NULL` (with a warning) when CLs cannot be
#' computed. Shared by the aggregate and decomposed se.fit paths.
#'
#' @param vcov_mat The fitted vcov (or NULL / wrong shape).
#' @param p Length of the parameter vector.
#' @return `list(vcov_use, free_idx)`, or `NULL` if unusable.
#' @keywords internal
.hzr_free_vcov <- function(vcov_mat, p) {
  vcov_ok <- !is.null(vcov_mat) && is.matrix(vcov_mat) &&
               nrow(vcov_mat) == p && ncol(vcov_mat) == p
  if (!vcov_ok) {
    warning("Variance-covariance matrix is unavailable; ",
            "standard errors and CLs will be NA.", call. = FALSE)
    return(NULL)
  }
  free_idx <- which(is.finite(diag(vcov_mat)))
  if (length(free_idx) < p) {
    free_submat <- vcov_mat[free_idx, free_idx, drop = FALSE]
    if (anyNA(free_submat)) {
      warning("Variance-covariance matrix has NA entries outside the ",
              "fixed-parameter rows/cols; standard errors and CLs will be NA.",
              call. = FALSE)
      return(NULL)
    }
    return(list(vcov_use = free_submat, free_idx = free_idx))
  }
  if (anyNA(vcov_mat)) {
    warning("Variance-covariance matrix has NA entries; ",
            "standard errors and CLs will be NA.", call. = FALSE)
    return(NULL)
  }
  list(vcov_use = vcov_mat, free_idx = free_idx)
}
```

Then in `.hzr_predict_with_se`, replace the whole vcov block — from the comment `# Fixed parameters (from \`fixed = "shapes"\` ...` through the `} else { vcov_use <- vcov_mat }` (the original lines that compute `vcov_ok`, `diag_vcov`, `free_idx`, `vcov_use`, and the three NA-return branches) — with:

```r
  fv <- .hzr_free_vcov(object$fit$vcov, p)
  if (is.null(fv)) {
    n <- length(target)
    fit <- if (type == "survival") exp(-target) else target
    na_vec <- rep(NA_real_, n)
    return(data.frame(fit = fit, se.fit = na_vec,
                      lower = na_vec, upper = na_vec))
  }
  vcov_use <- fv$vcov_use
  free_idx <- fv$free_idx
```

Leave the rest of `.hzr_predict_with_se` (Jacobian build, `J <- J[, free_idx, drop = FALSE]`, sandwich, scale transform) unchanged.

- [ ] **Step 4: Run tests; expect PASS**

Run: `Rscript -e 'suppressMessages(devtools::load_all(quiet=TRUE)); testthat::test_file("tests/testthat/test-predict-cl.R")'`
Expected: PASS — the new helper test AND every pre-existing se.fit test (notably the "vcov NA returns NA CLs with a warning" and "fixed shapes / CoE" tests) stay green, proving the refactor is behavior-preserving.

- [ ] **Step 5: Commit**

```bash
git add R/predict-cl.R tests/testthat/test-predict-cl.R
git commit -m "refactor(predict-cl): extract shared .hzr_free_vcov() helper

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: Decomposed se.fit driver `.hzr_predict_with_se_decomposed()`

**Files:**
- Modify: `R/predict-cl.R` (new function)
- Test: `tests/testthat/test-predict-cl.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-predict-cl.R`:

```r
test_that(".hzr_predict_with_se_decomposed returns a tidy long frame", {
  skip_on_cran()
  df <- make_toy()
  phases <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes"),
    constant = hzr_phase("constant")
  )
  fit <- hazard(survival::Surv(time, status) ~ 1, data = df,
                dist = "multiphase", phases = phases, fit = TRUE,
                control = list(n_starts = 2))
  t_new <- c(0.5, 1, 2)
  cov_counts <- c(early = 0L, constant = 0L)
  x_list <- list(early = NULL, constant = NULL)

  res <- TemporalHazard:::.hzr_predict_with_se_decomposed(
    object = fit, time = t_new, x_list = x_list,
    cov_counts = cov_counts, phases = phases, level = 0.95
  )

  expect_s3_class(res, "data.frame")
  expect_named(res, c("time", "component", "fit", "se.fit", "lower", "upper"))
  expect_equal(levels(res$component), c("total", "early", "constant"))
  expect_equal(nrow(res), length(t_new) * 3L)  # total + 2 phases
  # Monotone, non-negative CLs everywhere SEs are finite.
  ok <- is.finite(res$se.fit)
  expect_true(all(res$lower[ok] <= res$fit[ok] + 1e-9))
  expect_true(all(res$fit[ok] <= res$upper[ok] + 1e-9))
  expect_true(all(res$lower[ok] >= -1e-12))
})
```

- [ ] **Step 2: Run it; expect FAIL**

Run: `Rscript -e 'suppressMessages(devtools::load_all(quiet=TRUE)); testthat::test_file("tests/testthat/test-predict-cl.R")'`
Expected: FAIL — `.hzr_predict_with_se_decomposed` not found.

- [ ] **Step 3: Implement the driver**

Add to `R/predict-cl.R` (after `.hzr_predict_with_se`):

```r
#' Per-phase + total cumulative-hazard predictions with delta-method CLs
#'
#' Long-format companion to [.hzr_predict_with_se()] for multiphase models under
#' `decompose = TRUE`. Computes log-scale CLs for the total cumulative hazard
#' and for each phase's additive contribution `H_j(t) = mu_j(x) Phi_j(t)`.
#' Per-phase CLs use only that phase's parameter block, so they do NOT sum to
#' the total CL (cross-phase covariance contributes only to the total).
#'
#' @param object Fitted multiphase `hazard` object.
#' @param time Numeric prediction times (length n).
#' @param x_list Named list of per-phase design matrices.
#' @param cov_counts Named integer covariate counts.
#' @param phases Named list of `hzr_phase` objects.
#' @param level Confidence level.
#' @return Long `data.frame(time, component, fit, se.fit, lower, upper)`;
#'   `component` is an ordered factor with levels `c("total", names(phases))`.
#' @keywords internal
.hzr_predict_with_se_decomposed <- function(object, time, x_list,
                                             cov_counts, phases,
                                             level = 0.95) {
  theta <- object$fit$theta
  p <- length(theta)
  components <- c("total", names(phases))

  # Point estimates: total + per-phase cumulative hazard.
  ph <- .hzr_multiphase_cumhaz(time, theta, phases, cov_counts, x_list,
                               per_phase = TRUE)
  fit_by_comp <- c(list(total = ph$total),
                   stats::setNames(lapply(names(phases), function(nm) ph[[nm]]),
                                   names(phases)))

  make_long <- function(cl_list) {
    rows <- lapply(components, function(cmp) {
      data.frame(time = time, component = cmp,
                 fit = cl_list[[cmp]]$fit, se.fit = cl_list[[cmp]]$se.fit,
                 lower = cl_list[[cmp]]$lower, upper = cl_list[[cmp]]$upper,
                 stringsAsFactors = FALSE)
    })
    res <- do.call(rbind, rows)
    res$component <- factor(res$component, levels = components)
    rownames(res) <- NULL
    res
  }

  fv <- .hzr_free_vcov(object$fit$vcov, p)
  if (is.null(fv)) {
    na_cl <- lapply(components, function(cmp) {
      n <- length(time)
      data.frame(fit = fit_by_comp[[cmp]], se.fit = rep(NA_real_, n),
                 lower = rep(NA_real_, n), upper = rep(NA_real_, n))
    })
    names(na_cl) <- components
    return(make_long(na_cl))
  }

  # Per-phase Jacobians; the total is their sum.
  J_list <- .hzr_predict_jacobian_multiphase(theta, time, phases, cov_counts,
                                             x_list, p, per_phase = TRUE)
  J_by_comp <- c(list(total = Reduce(`+`, J_list)), J_list)

  cl_list <- lapply(components, function(cmp) {
    J <- J_by_comp[[cmp]][, fv$free_idx, drop = FALSE]
    se <- .hzr_predict_se_from_jacobian(J, fv$vcov_use)
    .hzr_predict_cl_from_se(fit_by_comp[[cmp]], se, level, "log")
  })
  names(cl_list) <- components
  make_long(cl_list)
}
```

- [ ] **Step 4: Run it; expect PASS**

Run: `Rscript -e 'suppressMessages(devtools::load_all(quiet=TRUE)); testthat::test_file("tests/testthat/test-predict-cl.R")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/predict-cl.R tests/testthat/test-predict-cl.R
git commit -m "feat(predict-cl): decomposed per-phase se.fit driver

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4: Wire into `predict.hazard()` + correctness tests

**Files:**
- Modify: `R/hazard_api.R` (the `if (decompose)` stop ~line 701; the multiphase `if (se.fit)` block ~line 861)
- Modify: `tests/testthat/test-predict-cl.R` (replace the line-58 error test)

- [ ] **Step 1: Write the end-to-end tests (and replace the stale error test)**

In `tests/testthat/test-predict-cl.R`, REPLACE the existing block:

```r
test_that("se.fit = TRUE with decompose = TRUE errors clearly", {
  df <- make_toy()
  phases <- list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                     fixed = "shapes"),
                 constant = hzr_phase("constant"))
  fit <- hazard(survival::Surv(time, status) ~ 1, data = df,
                dist = "multiphase", phases = phases, fit = TRUE,
                control = list(n_starts = 2))
  expect_error(
    predict(fit, newdata = data.frame(time = c(0.5, 1)),
            type = "cumulative_hazard", se.fit = TRUE, decompose = TRUE),
    "decompose"
  )
})
```

with:

```r
test_that("decompose + se.fit: cumulative_hazard works, survival errors", {
  skip_on_cran()
  df <- make_toy()
  phases <- list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                     fixed = "shapes"),
                 constant = hzr_phase("constant"))
  fit <- hazard(survival::Surv(time, status) ~ 1, data = df,
                dist = "multiphase", phases = phases, fit = TRUE,
                control = list(n_starts = 2))
  t_new <- c(0.5, 1, 2)

  # cumulative_hazard: long per-phase + total frame.
  res <- predict(fit, newdata = data.frame(time = t_new),
                 type = "cumulative_hazard", se.fit = TRUE, decompose = TRUE)
  expect_named(res, c("time", "component", "fit", "se.fit", "lower", "upper"))
  expect_equal(levels(res$component), c("total", "early", "constant"))
  expect_equal(nrow(res), length(t_new) * 3L)

  # Total rows reproduce the aggregate se.fit output exactly.
  agg <- predict(fit, newdata = data.frame(time = t_new),
                 type = "cumulative_hazard", se.fit = TRUE, decompose = FALSE)
  tot <- res[res$component == "total", ]
  expect_equal(tot$fit, agg$fit, tolerance = 1e-10)
  expect_equal(tot$se.fit, agg$se.fit, tolerance = 1e-10)
  expect_equal(tot$lower, agg$lower, tolerance = 1e-10)
  expect_equal(tot$upper, agg$upper, tolerance = 1e-10)

  # Per-phase point estimates reproduce the decompose=FALSE wide cumhaz columns.
  wide <- predict(fit, newdata = data.frame(time = t_new),
                  type = "cumulative_hazard", se.fit = FALSE, decompose = TRUE)
  expect_equal(res[res$component == "early", "fit"], wide$early, tolerance = 1e-10)
  expect_equal(res[res$component == "constant", "fit"], wide$constant,
               tolerance = 1e-10)

  # survival still rejected.
  expect_error(
    predict(fit, newdata = data.frame(time = t_new),
            type = "survival", se.fit = TRUE, decompose = TRUE),
    "cumulative_hazard"
  )
})

test_that("decompose + se.fit per-phase SE matches a numeric jacobian", {
  skip_if_not_installed("numDeriv")
  skip_on_cran()
  df <- make_toy()
  phases <- list(early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1,
                                     fixed = "shapes"),
                 constant = hzr_phase("constant"))
  fit <- hazard(survival::Surv(time, status) ~ 1, data = df,
                dist = "multiphase", phases = phases, fit = TRUE,
                control = list(n_starts = 2))
  t_new <- c(0.5, 1, 2)
  cov_counts <- c(early = 0L, constant = 0L)
  x_list <- list(early = NULL, constant = NULL)
  theta <- fit$fit$theta
  V <- fit$fit$vcov
  free_idx <- which(is.finite(diag(V)))

  res <- predict(fit, newdata = data.frame(time = t_new),
                 type = "cumulative_hazard", se.fit = TRUE, decompose = TRUE)

  # Numeric per-phase SE for the "early" component via numDeriv.
  H_early <- function(th) {
    TemporalHazard:::.hzr_multiphase_cumhaz(
      t_new, th, phases, cov_counts, x_list, per_phase = TRUE
    )$early
  }
  Jn <- numDeriv::jacobian(H_early, theta)[, free_idx, drop = FALSE]
  se_num <- sqrt(diag(Jn %*% V[free_idx, free_idx, drop = FALSE] %*% t(Jn)))
  se_ana <- res[res$component == "early", "se.fit"]
  expect_equal(se_ana, se_num, tolerance = 1e-5)
})
```

- [ ] **Step 2: Run them; expect FAIL**

Run: `Rscript -e 'suppressMessages(devtools::load_all(quiet=TRUE)); testthat::test_file("tests/testthat/test-predict-cl.R")'`
Expected: FAIL — `predict(... decompose=TRUE, se.fit=TRUE, type="cumulative_hazard")` still errors with the old "decompose" message.

- [ ] **Step 3: Update the validation stop**

In `R/hazard_api.R`, replace the existing block:

```r
    if (decompose) {
      stop("'se.fit = TRUE' with 'decompose = TRUE' is not supported. ",
           "Request point predictions first (decompose = TRUE, ",
           "se.fit = FALSE), then compute CLs separately on the total.",
           call. = FALSE)
    }
```

with:

```r
    if (decompose && type != "cumulative_hazard") {
      stop("'se.fit = TRUE' with 'decompose = TRUE' is only supported for ",
           "type = \"cumulative_hazard\" (per-phase survival is not additive). ",
           "Got type = \"", type, "\".", call. = FALSE)
    }
```

- [ ] **Step 4: Dispatch in the multiphase branch**

In `R/hazard_api.R`, inside the multiphase prediction branch, replace:

```r
      if (se.fit) {
        diff_fn <- function(th) {
          .hzr_multiphase_cumhaz(pred_time, th, phases, cov_counts, x_list)
        }
        return(.hzr_predict_with_se(
          object = object, type = type, time = pred_time,
          x_list = x_list, cov_counts = cov_counts, phases = phases,
          level = level, diff_fn = diff_fn
        ))
      }
```

with:

```r
      if (se.fit) {
        if (decompose) {
          return(.hzr_predict_with_se_decomposed(
            object = object, time = pred_time, x_list = x_list,
            cov_counts = cov_counts, phases = phases, level = level
          ))
        }
        diff_fn <- function(th) {
          .hzr_multiphase_cumhaz(pred_time, th, phases, cov_counts, x_list)
        }
        return(.hzr_predict_with_se(
          object = object, type = type, time = pred_time,
          x_list = x_list, cov_counts = cov_counts, phases = phases,
          level = level, diff_fn = diff_fn
        ))
      }
```

- [ ] **Step 5: Run tests; expect PASS**

Run: `Rscript -e 'suppressMessages(devtools::load_all(quiet=TRUE)); testthat::test_file("tests/testthat/test-predict-cl.R")'`
Expected: PASS (all new + pre-existing tests).

- [ ] **Step 6: Commit**

```bash
git add R/hazard_api.R tests/testthat/test-predict-cl.R
git commit -m "feat(predict): per-phase CLs for decompose=TRUE, se.fit=TRUE (cumulative_hazard)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 5: Docs, NEWS, and full verification

**Files:**
- Modify: `R/hazard_api.R` (roxygen `@param decompose` / `@param se.fit` / `@details`, ~lines 528-547)
- Modify: `NEWS.md`
- Regenerate: `man/predict.hazard.Rd` (via `devtools::document()`)

- [ ] **Step 1: Update roxygen on `predict.hazard()`**

In `R/hazard_api.R`, find the `@param se.fit` / `@param decompose` / `@details` text that says the combination is unsupported (the lines stating "Not compatible with `decompose = TRUE`" and "`decompose = TRUE` and `se.fit = TRUE` cannot be used together ..."). Replace that guidance with:

```r
#'   For multiphase models, `se.fit = TRUE` combines with `decompose = TRUE`
#'   when `type = "cumulative_hazard"`: the result is a long data frame with one
#'   row per prediction time and component (`component` in `"total"` plus each
#'   phase name) and columns `fit`, `se.fit`, `lower`, `upper`. Per-phase CLs
#'   use only that phase's parameters, so they do not sum to the total CL.
#'   The combination is not available for `type = "survival"` (per-phase
#'   survival is not additive).
```

(Keep the `@param`/`@details` structure; only swap the now-inaccurate "cannot be used together" prose.)

- [ ] **Step 2: Add a NEWS entry**

In `NEWS.md`, under the top `## New features` heading (create the heading just under the version header if it does not exist, above `## Bug fixes`), add:

```md
* `predict.hazard(type = "cumulative_hazard", decompose = TRUE, se.fit = TRUE)`
  now returns per-phase **and** total delta-method confidence limits for
  multiphase models, as a long data frame
  (`time`, `component`, `fit`, `se.fit`, `lower`, `upper`). Each phase's CL uses
  only that phase's parameters, so per-phase limits do not sum to the total.
  Previously this combination raised an error.
```

- [ ] **Step 3: Regenerate Rd and run full verification**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'devtools::test()'
Rscript -e 'lintr::lint("R/predict-cl.R"); lintr::lint("R/hazard_api.R")'
```
Expected: `document()` updates `man/predict.hazard.Rd` with no errors; `devtools::test()` reports 0 failures; `lintr` reports no lints.

- [ ] **Step 4: Commit**

```bash
git add R/hazard_api.R NEWS.md man/predict.hazard.Rd
git commit -m "docs: document per-phase CLs for predict(decompose, se.fit)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

- [ ] **Step 5: Push and open PR**

```bash
git push -u origin feat/predict-decompose-se
gh pr create --base dev --title "feat(predict): per-phase CLs for decompose=TRUE, se.fit=TRUE (7c)" \
  --body "Unblocks predict(type='cumulative_hazard', decompose=TRUE, se.fit=TRUE) for multiphase models. Returns long per-phase + total delta-method CLs. Reuses the existing per-phase Jacobian blocks + delta-method sandwich (shared via new .hzr_free_vcov() helper). Spec + plan in inst/dev/."
```

---

## Notes for the implementer

- If `devtools` is unavailable, substitute `R CMD INSTALL . && Rscript -e 'testthat::test_dir("tests/testthat")'` and `roxygen2::roxygenise()` for `document()`.
- Tolerances mirror the surrounding tests: `1e-10` for exact reproductions (total == aggregate, per-phase == wide), `1e-5` for analytic-vs-numeric SE, `1e-6`/`1e-12` for Jacobian checks.
- Do not loosen a tolerance or relax a shape assertion to make a test pass — a failure is a real defect.
- The per-phase-vs-numeric SE test (Task 4) is the core correctness anchor for the new Jacobian-block math; if it fails, debug the `per_phase` column placement in Task 1, not the test.
