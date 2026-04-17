# stepwise-step.R — Per-step drivers for hzr_stepwise()
#
# Step 8.4–8.5 of STEPWISE-DESIGN.md.  One forward step or one backward
# step, driven by the unified `.hzr_candidate_score()` wrapper.  Higher
# level control (two-way loop, MOVE cap, force_in, max_steps) lives in
# the hzr_stepwise() driver to be built in §8.6.

# ---------------------------------------------------------------------------
# Candidate enumeration
# ---------------------------------------------------------------------------

#' Normalise a user-supplied scope into a flat list of (var, phase) pairs
#'
#' Drops variables already in the model, variables in `force_out`, and
#' candidates whose phase is not present in a multiphase fit.
#'
#' @param fit The current fit.
#' @param scope One of: NULL (data-driven; infer from `data`), a
#'   one-sided formula (single-dist), or a named list of one-sided
#'   formulas keyed by phase (multiphase).
#' @param data Data frame to derive candidates from when `scope = NULL`.
#' @param force_out Character vector of variable names to exclude.
#'
#' @return A list of `list(var = <chr>, phase = <chr or NULL>)`.
#'
#' @keywords internal
#' @noRd
.hzr_stepwise_candidates <- function(fit, scope = NULL, data,
                                      force_out = character()) {
  dist <- fit$spec$dist

  if (dist == "multiphase") {
    phase_names <- names(fit$spec$phases)
    current_per_phase <- .hzr_scope_current_vars(fit)

    if (is.null(scope)) {
      # Default: all data-frame vars (excluding Surv components and
      # force_out) are candidates for every phase.
      lhs_vars <- all.vars(stats::as.formula(
        deparse(fit$call$formula)
      )[[2L]])
      data_vars <- setdiff(colnames(data), c(lhs_vars, force_out))
      scope <- setNames(
        lapply(phase_names, function(p) {
          rhs_syms <- data_vars
          # Return as a one-sided formula for symmetry with the
          # user-supplied case.
          if (length(rhs_syms) == 0L) return(NULL)
          stats::as.formula(paste("~", paste(rhs_syms, collapse = " + ")))
        }),
        phase_names
      )
    } else {
      if (!is.list(scope) || is.null(names(scope))) {
        stop("`scope` for multiphase fits must be a named list keyed by phase.",
             call. = FALSE)
      }
      unknown <- setdiff(names(scope), phase_names)
      if (length(unknown) > 0L) {
        stop("Unknown phase(s) in scope: ",
             paste(sQuote(unknown), collapse = ", "),
             call. = FALSE)
      }
    }

    candidates <- list()
    for (p in names(scope)) {
      sc <- scope[[p]]
      if (is.null(sc)) next
      terms_p <- .hzr_formula_rhs_terms(sc)
      eligible <- setdiff(terms_p, c(current_per_phase[[p]], force_out))
      for (v in eligible) {
        candidates[[length(candidates) + 1L]] <-
          list(var = v, phase = p)
      }
    }
    return(candidates)
  }

  # Single-distribution
  if (is.null(scope)) {
    f <- fit$call$formula
    if (!is.null(f) && !inherits(f, "formula")) {
      f <- stats::as.formula(deparse(f))
    }
    lhs_vars <- if (is.null(f)) character() else all.vars(f[[2L]])
    data_vars <- setdiff(colnames(data), c(lhs_vars, force_out))
  } else {
    if (inherits(scope, "formula")) {
      data_vars <- .hzr_formula_rhs_terms(scope)
    } else if (is.character(scope)) {
      data_vars <- scope
    } else {
      stop("`scope` must be NULL, a one-sided formula, or a character vector.",
           call. = FALSE)
    }
    data_vars <- setdiff(data_vars, force_out)
  }

  current_vars <- .hzr_scope_current_vars(fit)
  eligible <- setdiff(data_vars, current_vars)
  lapply(eligible, function(v) list(var = v, phase = NULL))
}


# ---------------------------------------------------------------------------
# Forward step
# ---------------------------------------------------------------------------

#' Execute one forward step of stepwise selection
#'
#' Refits the current model with each candidate (variable, phase)
#' appended, scores each via `.hzr_candidate_score(mode = "entry")`,
#' and accepts the best if its score clears the entry threshold.
#'
#' Divergent candidate refits emit a `warning()` naming the failing
#' `(variable, phase)` pair and are excluded from the selection, per
#' the §2 Q5 decision in STEPWISE-DESIGN.md.
#'
#' @param current Fitted `hazard` object that is the starting point.
#' @param scope Scope specification (see `.hzr_stepwise_candidates`).
#' @param data Data frame for refits.
#' @param criterion Either `"wald"` or `"aic"`.
#' @param slentry Entry threshold for the Wald criterion (ignored when
#'   `criterion = "aic"`; the entry rule there is ΔAIC < 0).
#' @param force_out Character vector of variables that may never be
#'   considered as candidates.
#' @param ... Forwarded to `.hzr_refit_with_scope()` (and thence to
#'   `hazard()`), e.g. `control = list(...)`.
#'
#' @return A list with:
#' \describe{
#'   \item{accepted}{Logical; did any candidate enter this step?}
#'   \item{fit}{If accepted, the new fit.  Otherwise `current` echoed.}
#'   \item{variable}{Variable that entered, or `NA_character_`.}
#'   \item{phase}{Phase entered, or `NA_character_` (single-dist).}
#'   \item{score}{Winning score (p or ΔAIC), or `NA_real_`.}
#'   \item{p_value}{Winning p-value.}
#'   \item{delta_aic}{Winning ΔAIC.}
#'   \item{stat, df}{Wald statistic / df of the winner.}
#'   \item{all_scores}{Tibble-like data frame of every candidate
#'     considered and its score.}
#'   \item{refit_failures}{Character vector of `"var@phase"` tokens for
#'     candidates whose refit diverged.}
#' }
#'
#' @keywords internal
#' @noRd
.hzr_stepwise_forward_step <- function(current, scope = NULL, data,
                                        criterion = c("wald", "aic"),
                                        slentry   = 0.30,
                                        force_out = character(),
                                        ...) {
  criterion <- match.arg(criterion)
  if (!inherits(current, "hazard")) {
    stop("`current` must be a fitted `hazard` object.", call. = FALSE)
  }

  cands <- .hzr_stepwise_candidates(current, scope = scope, data = data,
                                     force_out = force_out)

  null_result <- function() {
    list(
      accepted  = FALSE,
      fit       = current,
      variable  = NA_character_,
      phase     = NA_character_,
      score     = NA_real_,
      p_value   = NA_real_,
      delta_aic = NA_real_,
      stat      = NA_real_,
      df        = NA_integer_,
      all_scores = data.frame(
        variable  = character(),
        phase     = character(),
        score     = numeric(),
        p_value   = numeric(),
        delta_aic = numeric(),
        stat      = numeric(),
        df        = integer(),
        stringsAsFactors = FALSE
      ),
      refit_failures = character()
    )
  }

  if (length(cands) == 0L) {
    return(null_result())
  }

  rows     <- vector("list", length(cands))
  failures <- character()

  for (i in seq_along(cands)) {
    cand <- cands[[i]]
    candidate_fit <- tryCatch(
      .hzr_refit_with_scope(
        current, action = "add",
        var = cand$var, phase = cand$phase,
        data = data, ...
      ),
      error = function(e) e
    )

    failure_token <- if (is.null(cand$phase)) {
      cand$var
    } else {
      paste0(cand$var, "@", cand$phase)
    }

    if (inherits(candidate_fit, "error") ||
          isFALSE(candidate_fit$fit$converged)) {
      warning("Stepwise forward: candidate refit failed for ",
              failure_token, ".", call. = FALSE)
      failures <- c(failures, failure_token)
      rows[[i]] <- data.frame(
        variable  = cand$var,
        phase     = cand$phase %||% NA_character_,
        score     = NA_real_,
        p_value   = NA_real_,
        delta_aic = NA_real_,
        stat      = NA_real_,
        df        = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }

    # Coefficient name(s) of the newly-entered variable in the candidate
    # fit.  For single-dist this is `var` directly; for multiphase the
    # coef is phase-prefixed.
    coef_name <- .hzr_candidate_coef_name(candidate_fit, cand$var,
                                           cand$phase)

    s <- .hzr_candidate_score(
      criterion = criterion, mode = "entry",
      current = current, candidate = candidate_fit,
      names = coef_name
    )

    rows[[i]] <- data.frame(
      variable  = cand$var,
      phase     = cand$phase %||% NA_character_,
      score     = s$score,
      p_value   = s$p_value,
      delta_aic = s$delta_aic,
      stat      = s$stat,
      df        = s$df,
      stringsAsFactors = FALSE
    )
    # Cache the fit on the row so we can recover it for the winner
    attr(rows[[i]], "fit") <- candidate_fit
  }

  all_scores <- do.call(rbind, rows)
  # Strip the per-row fit attributes from the combined frame but keep
  # them in a parallel list keyed by row for winner lookup.
  candidate_fits <- lapply(rows, function(r) attr(r, "fit"))

  valid <- which(!is.na(all_scores$score))
  if (length(valid) == 0L) {
    out <- null_result()
    out$all_scores <- all_scores
    out$refit_failures <- failures
    return(out)
  }

  best_idx <- valid[which.min(all_scores$score[valid])]
  best     <- all_scores[best_idx, ]

  threshold_met <- if (criterion == "wald") {
    best$score < slentry
  } else {
    best$score < 0
  }

  if (!threshold_met) {
    out <- null_result()
    out$all_scores <- all_scores
    out$refit_failures <- failures
    return(out)
  }

  list(
    accepted  = TRUE,
    fit       = candidate_fits[[best_idx]],
    variable  = best$variable,
    phase     = best$phase,
    score     = best$score,
    p_value   = best$p_value,
    delta_aic = best$delta_aic,
    stat      = best$stat,
    df        = best$df,
    all_scores = all_scores,
    refit_failures = failures
  )
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Backward step
# ---------------------------------------------------------------------------

#' Enumerate the currently-in-model (variable, phase) pairs eligible
#' for a drop test.
#'
#' Force-in variables are *included* in the returned list so their
#' scores appear in the trace (§2 Q4 of STEPWISE-DESIGN.md); the
#' forward/backward driver filters them out of the argmax pool via the
#' `force_in` flag on each row.
#'
#' @keywords internal
#' @noRd
.hzr_stepwise_drop_candidates <- function(fit) {
  dist <- fit$spec$dist
  if (dist == "multiphase") {
    per_phase <- .hzr_scope_current_vars(fit)
    out <- list()
    for (p in names(per_phase)) {
      for (v in per_phase[[p]]) {
        out[[length(out) + 1L]] <- list(var = v, phase = p)
      }
    }
    return(out)
  }
  vars <- .hzr_scope_current_vars(fit)
  lapply(vars, function(v) list(var = v, phase = NULL))
}


#' Execute one backward step of stepwise selection
#'
#' For each currently-in-model variable, scores drop via
#' `.hzr_candidate_score(mode = "drop")`.  Force-in variables are
#' scored and reported but excluded from the drop decision.
#'
#' Unlike the forward step there is no per-candidate refit: the Wald
#' statistic and its AIC approximation both operate on the current
#' model's vcov.  A single refit fires only after the drop decision
#' via `.hzr_refit_with_scope()`.
#'
#' @param current Fitted `hazard` object.
#' @param data Data frame used to rebuild the dropped model.
#' @param criterion Either `"wald"` or `"aic"`.
#' @param slstay Retention threshold for the Wald criterion (ignored
#'   when `criterion = "aic"`; the drop rule there is ΔAIC_drop < 0).
#' @param force_in Character vector of variables that may never be
#'   dropped.
#' @param ... Forwarded to `.hzr_refit_with_scope()` for the post-drop
#'   refit.
#'
#' @return A list with the same top-level shape as
#'   `.hzr_stepwise_forward_step()`, except:
#'   * `all_scores` gains a logical `force_in` column.
#'   * The action this represents is a drop, so `accepted = TRUE` means
#'     the variable was removed from the model.
#'
#' @keywords internal
#' @noRd
.hzr_stepwise_backward_step <- function(current, data,
                                         criterion = c("wald", "aic"),
                                         slstay    = 0.20,
                                         force_in  = character(),
                                         ...) {
  criterion <- match.arg(criterion)
  if (!inherits(current, "hazard")) {
    stop("`current` must be a fitted `hazard` object.", call. = FALSE)
  }

  cands <- .hzr_stepwise_drop_candidates(current)

  empty_scores <- data.frame(
    variable  = character(),
    phase     = character(),
    force_in  = logical(),
    score     = numeric(),
    p_value   = numeric(),
    delta_aic = numeric(),
    stat      = numeric(),
    df        = integer(),
    stringsAsFactors = FALSE
  )

  null_result <- function(all_scores = empty_scores) {
    list(
      accepted  = FALSE,
      fit       = current,
      variable  = NA_character_,
      phase     = NA_character_,
      score     = NA_real_,
      p_value   = NA_real_,
      delta_aic = NA_real_,
      stat      = NA_real_,
      df        = NA_integer_,
      all_scores     = all_scores,
      refit_failures = character()
    )
  }

  if (length(cands) == 0L) {
    return(null_result())
  }

  rows <- vector("list", length(cands))
  for (i in seq_along(cands)) {
    cand <- cands[[i]]
    coef_name <- .hzr_candidate_coef_name(current, cand$var, cand$phase)

    s <- .hzr_candidate_score(
      criterion = criterion, mode = "drop",
      current = current, names = coef_name
    )

    rows[[i]] <- data.frame(
      variable  = cand$var,
      phase     = cand$phase %||% NA_character_,
      force_in  = cand$var %in% force_in,
      score     = s$score,
      p_value   = s$p_value,
      delta_aic = s$delta_aic,
      stat      = s$stat,
      df        = s$df,
      stringsAsFactors = FALSE
    )
  }
  all_scores <- do.call(rbind, rows)

  eligible <- which(!all_scores$force_in & !is.na(all_scores$score))
  if (length(eligible) == 0L) {
    return(null_result(all_scores))
  }

  best_idx <- eligible[which.min(all_scores$score[eligible])]
  best     <- all_scores[best_idx, ]

  threshold_met <- if (criterion == "wald") {
    best$score < (1 - slstay)          # i.e. p > slstay
  } else {
    best$score < 0                     # ΔAIC_drop < 0
  }

  if (!threshold_met) {
    return(null_result(all_scores))
  }

  refitted <- tryCatch(
    .hzr_refit_with_scope(
      current, action = "drop",
      var = best$variable,
      phase = if (is.na(best$phase)) NULL else best$phase,
      data = data, ...
    ),
    error = function(e) e
  )

  failure_token <- if (is.na(best$phase)) {
    best$variable
  } else {
    paste0(best$variable, "@", best$phase)
  }

  if (inherits(refitted, "error") || isFALSE(refitted$fit$converged)) {
    warning("Stepwise backward: post-drop refit failed for ",
            failure_token, ".", call. = FALSE)
    out <- null_result(all_scores)
    out$refit_failures <- failure_token
    return(out)
  }

  list(
    accepted  = TRUE,
    fit       = refitted,
    variable  = best$variable,
    phase     = best$phase,
    score     = best$score,
    p_value   = best$p_value,
    delta_aic = best$delta_aic,
    stat      = best$stat,
    df        = best$df,
    all_scores     = all_scores,
    refit_failures = character()
  )
}


#' Name under which a newly-entered variable appears in coef(fit)
#'
#' Canonical naming differs between fit kinds:
#'   multiphase  — phase-prefixed formula names (e.g. `"early.age"`).
#'   single-dist — positional `"betaN"` from `.hzr_parameter_names()`,
#'     where N is the column index of `var` in `colnames(fit$data$x)`.
#'
#' This matches the naming `summary.hazard()` prints and the canonical
#' name `.hzr_wald_p()` uses for coefficient lookup.
#'
#' @keywords internal
#' @noRd
.hzr_candidate_coef_name <- function(fit, var, phase) {
  if (fit$spec$dist == "multiphase") {
    return(paste0(phase, ".", var))
  }
  xcols <- colnames(fit$data$x)
  idx <- match(var, xcols)
  if (is.na(idx)) {
    stop("Variable ", sQuote(var),
         " not found in the design matrix.",
         call. = FALSE)
  }
  paste0("beta", idx)
}


#' %||% — NULL-coalesce for lazy defaults
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
