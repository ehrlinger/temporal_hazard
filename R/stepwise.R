# stepwise.R — user-facing driver for Phase 4b stepwise covariate
# selection.  Step 8.6 of STEPWISE-DESIGN.md.
#
# The driver combines .hzr_stepwise_forward_step() and
# .hzr_stepwise_backward_step() into a two-way loop with:
#
#   * SAS-style MOVE oscillation guard: when a variable hits
#     `max_move` entries/exits, it is frozen (added to both the
#     force_in and force_out sets used internally).
#   * `max_steps` hard cap that warns on hit.
#   * `force_in` / `force_out` user-supplied constraints.
#   * Optional console trace mirroring the §5 output spec.
#
# The returned object inherits from both `hzr_stepwise` and `hazard`
# so it can be passed to `predict()`, `summary()`, `coef()`, etc.
# through the same S3 infrastructure as a plain fit.

#' Stepwise covariate selection for a parametric hazard model
#'
#' Run forward, backward, or two-way stepwise selection on an existing
#' `hazard` fit using Wald p-values or AIC deltas as the entry /
#' retention criterion.  Phase-specific entry is supported for
#' multiphase models: a covariate can enter one phase and not another.
#'
#' @param fit A fitted `hazard` object built via the
#'   `formula = Surv(...) ~ predictors, data = df` interface.
#' @param scope Candidate set.  `NULL` (default) uses every data-frame
#'   column not already in the model for every phase.  For
#'   single-distribution fits, pass a one-sided formula
#'   (`~ age + nyha`) or a character vector of names.  For multiphase
#'   fits, pass a named list of one-sided formulas keyed by phase.
#' @param data Data frame the base fit was built on.  Required for
#'   refits.
#' @param direction One of `"both"` (default), `"forward"`,
#'   `"backward"`.
#' @param criterion One of `"wald"` (default) or `"aic"`.  SAS-style
#'   p-value thresholds apply to Wald; AIC uses `ΔAIC < 0` uniformly.
#' @param slentry Entry p-value threshold for the Wald criterion.
#'   Default `0.30` matches SAS `SLENTRY`.
#' @param slstay Retention p-value threshold for the Wald criterion.
#'   Default `0.20` matches SAS `SLSTAY`.
#' @param max_steps Hard cap on total accepted actions.  Emits a
#'   `warning()` if hit.  Default `50`.
#' @param max_move Per-variable oscillation cap.  When a variable has
#'   entered + exited more than `max_move` times it is frozen for the
#'   remainder of the run.  Default `4`.
#' @param force_in Character vector of variables that must remain in
#'   the model.  Such variables are still scored and reported in the
#'   selection trace, but are never dropped.
#' @param force_out Character vector of variables that may never be
#'   considered as candidates.
#' @param trace Logical; print step-by-step progress to the console.
#'   Default `TRUE`.
#' @param ... Passed to the underlying `hazard()` refits (e.g.
#'   `control = list(n_starts = 3)`).
#'
#' @return An object of class `c("hzr_stepwise", "hazard")` — the
#'   final fit augmented with:
#'   \describe{
#'     \item{\code{steps}}{Data frame with one row per accepted /
#'       frozen action; see Details.}
#'     \item{\code{scope}}{Record of the candidate scope, plus
#'       `force_in`, `force_out`, and the frozen set.}
#'     \item{\code{criteria}}{Named list of the threshold / direction
#'       settings actually applied.}
#'     \item{\code{trace_msg}}{Character vector of the trace lines,
#'       captured regardless of the `trace` flag.}
#'     \item{\code{elapsed}}{`difftime` from start to finish.}
#'     \item{\code{final_call}}{The call that produced this result.}
#'   }
#'
#' @details
#' The `steps` data frame has columns:
#'
#' \describe{
#'   \item{\code{step_num}}{Integer sequence starting at 1.}
#'   \item{\code{action}}{`"enter"`, `"drop"`, or `"frozen"`.}
#'   \item{\code{variable}}{Variable affected.}
#'   \item{\code{phase}}{Phase name (multiphase) or `NA_character_`.}
#'   \item{\code{criterion}}{`"wald"` or `"aic"`.}
#'   \item{\code{score}}{Winning score used for the decision.}
#'   \item{\code{stat}, \code{df}}{Wald statistic and degrees of
#'     freedom.}
#'   \item{\code{p_value}, \code{delta_aic}}{Always populated when
#'     computable, regardless of the active criterion.}
#'   \item{\code{logLik}, \code{aic}, \code{n_coef}}{Goodness-of-fit
#'     diagnostics of the model *after* this step.}
#' }
#'
#' @export
hzr_stepwise <- function(fit,
                         scope     = NULL,
                         data,
                         direction = c("both", "forward", "backward"),
                         criterion = c("wald", "aic"),
                         slentry   = 0.30,
                         slstay    = 0.20,
                         max_steps = 50L,
                         max_move  = 4L,
                         force_in  = character(),
                         force_out = character(),
                         trace     = TRUE,
                         ...) {
  direction <- match.arg(direction)
  criterion <- match.arg(criterion)

  if (!inherits(fit, "hazard")) {
    stop("`fit` must be a `hazard` object.", call. = FALSE)
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("`data` must be a data frame (typically the frame used for the base fit).",
         call. = FALSE)
  }

  ts_start <- Sys.time()
  call <- match.call()

  extra_args <- list(...)

  steps     <- list()
  trace_msg <- character()

  emit <- function(msg) {
    trace_msg[[length(trace_msg) + 1L]] <<- msg
    if (isTRUE(trace)) cat(msg, "\n", sep = "")
  }

  # Header line (mirrors design §5).
  header <- if (criterion == "wald") {
    sprintf(
      "Stepwise selection (direction = %s, criterion = wald, slentry = %.2f, slstay = %.2f)",
      direction, slentry, slstay
    )
  } else {
    sprintf(
      "Stepwise selection (direction = %s, criterion = aic)",
      direction
    )
  }
  emit(header)
  emit("")

  # Move counter: per-variable tally of entries + exits.  Use a named
  # list rather than a named integer vector — `lst[[missing]]` returns
  # NULL, whereas `vec[[missing]]` errors with "subscript out of bounds".
  move_counts <- list()
  frozen      <- character()

  current <- fit
  step_no <- 0L
  stopped_by_max_steps <- FALSE

  record_step <- function(action, out) {
    step_no <<- step_no + 1L
    row <- data.frame(
      step_num  = step_no,
      action    = action,
      variable  = out$variable,
      phase     = out$phase,
      criterion = criterion,
      score     = out$score,
      stat      = out$stat,
      df        = out$df,
      p_value   = out$p_value,
      delta_aic = out$delta_aic,
      logLik    = current$fit$objective   %||% NA_real_,
      aic       = .hzr_aic(current),
      n_coef    = length(current$fit$theta),
      stringsAsFactors = FALSE
    )
    steps[[length(steps) + 1L]] <<- row

    score_fmt <- if (criterion == "wald") {
      sprintf("p = %.3f", out$p_value)
    } else {
      sprintf("\u0394AIC = %+.2f", out$delta_aic)
    }
    phase_txt <- if (is.na(out$phase)) {
      ""
    } else if (action == "enter") {
      paste0("  into  ", out$phase)
    } else {
      paste0("  from  ", out$phase)
    }
    emit(sprintf(
      "Step %d: %-6s %s%s   (%s)",
      step_no, toupper(action), out$variable, phase_txt, score_fmt
    ))
  }

  record_freeze <- function(var, phase_hint = NA_character_) {
    step_no <<- step_no + 1L
    row <- data.frame(
      step_num  = step_no,
      action    = "frozen",
      variable  = var,
      phase     = phase_hint,
      criterion = criterion,
      score     = NA_real_,
      stat      = NA_real_,
      df        = NA_integer_,
      p_value   = NA_real_,
      delta_aic = NA_real_,
      logLik    = current$fit$objective   %||% NA_real_,
      aic       = .hzr_aic(current),
      n_coef    = length(current$fit$theta),
      stringsAsFactors = FALSE
    )
    steps[[length(steps) + 1L]] <<- row
    emit(sprintf(
      "Step %d: FROZEN %s   (exceeded max_move = %d; OSCILLATING)",
      step_no, var, max_move
    ))
  }

  bump_move <- function(var) {
    move_counts[[var]] <<- (move_counts[[var]] %||% 0L) + 1L
    if (move_counts[[var]] > max_move && !var %in% frozen) {
      frozen <<- c(frozen, var)
      record_freeze(var)
    }
  }

  # Main loop
  repeat {
    if (step_no >= max_steps) {
      warning("Stepwise selection hit max_steps = ", max_steps,
              "; stopping early.", call. = FALSE)
      stopped_by_max_steps <- TRUE
      break
    }

    add_happened  <- FALSE
    drop_happened <- FALSE

    effective_force_out <- unique(c(force_out, frozen))
    effective_force_in  <- unique(c(force_in,  frozen))

    if (direction %in% c("forward", "both")) {
      fwd <- do.call(.hzr_stepwise_forward_step, c(list(
        current   = current,
        scope     = scope,
        data      = data,
        criterion = criterion,
        slentry   = slentry,
        force_out = effective_force_out
      ), extra_args))

      if (fwd$accepted) {
        current <- fwd$fit
        record_step("enter", fwd)
        bump_move(fwd$variable)
        add_happened <- TRUE
      }
    }

    if (direction %in% c("backward", "both")) {
      bwd <- do.call(.hzr_stepwise_backward_step, c(list(
        current   = current,
        data      = data,
        criterion = criterion,
        slstay    = slstay,
        force_in  = effective_force_in
      ), extra_args))

      if (bwd$accepted) {
        current <- bwd$fit
        record_step("drop", bwd)
        bump_move(bwd$variable)
        drop_happened <- TRUE
      }
    }

    if (!add_happened && !drop_happened) {
      emit(sprintf("(no further action after %d step%s)",
                   step_no, if (step_no == 1L) "" else "s"))
      break
    }
  }

  elapsed <- difftime(Sys.time(), ts_start, units = "secs")

  emit("")
  emit(sprintf("Final model: %d covariate%s, logLik = %.2f, AIC = %.2f",
               max(0L, length(current$fit$theta) -
                     .hzr_stepwise_shape_count(current)),
               if (length(current$fit$theta) -
                     .hzr_stepwise_shape_count(current) == 1L) "" else "s",
               current$fit$objective %||% NA_real_,
               .hzr_aic(current)))

  steps_df <- if (length(steps) == 0L) {
    data.frame(
      step_num = integer(), action = character(),
      variable = character(), phase = character(),
      criterion = character(), score = numeric(),
      stat = numeric(), df = integer(),
      p_value = numeric(), delta_aic = numeric(),
      logLik = numeric(), aic = numeric(), n_coef = integer(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, steps)
  }

  result <- current
  result$steps      <- steps_df
  result$scope      <- list(
    candidates = scope,
    force_in   = force_in,
    force_out  = force_out,
    frozen     = frozen
  )
  result$criteria   <- list(
    direction = direction,
    criterion = criterion,
    slentry   = slentry,
    slstay    = slstay,
    max_steps = max_steps,
    max_move  = max_move,
    hit_max_steps = stopped_by_max_steps
  )
  result$trace_msg  <- trace_msg
  result$elapsed    <- elapsed
  result$final_call <- call

  class(result) <- unique(c("hzr_stepwise", class(result)))
  result
}


#' Print method for `hzr_stepwise`
#'
#' Shows the full selection trace and a brief final-model summary.
#'
#' @param x An `hzr_stepwise` object.
#' @param ... Unused.
#' @export
print.hzr_stepwise <- function(x, ...) {
  cat(paste(x$trace_msg, collapse = "\n"), "\n", sep = "")
  invisible(x)
}


#' Summary method for `hzr_stepwise`
#'
#' Produces the standard `summary.hazard()` on the final fit, with the
#' selection trace prepended in the print method.
#'
#' @param object An `hzr_stepwise` object.
#' @param ... Unused.
#' @export
summary.hzr_stepwise <- function(object, ...) {
  # Strip the stepwise class so NextMethod dispatches cleanly to
  # summary.hazard.
  class(object) <- setdiff(class(object), "hzr_stepwise")
  out <- NextMethod()
  out$stepwise_steps <- object$steps
  out$stepwise_trace <- object$trace_msg
  class(out) <- unique(c("summary.hzr_stepwise", class(out)))
  out
}


#' @export
print.summary.hzr_stepwise <- function(x, ...) {
  if (!is.null(x$stepwise_trace)) {
    cat(paste(x$stepwise_trace, collapse = "\n"), "\n\n", sep = "")
  }
  class(x) <- setdiff(class(x), "summary.hzr_stepwise")
  NextMethod()
}


#' Coerce an `hzr_stepwise` result to its selection trace
#'
#' Returns the `$steps` data frame so downstream tidyverse / data.table
#' pipelines can work with the trace directly.
#'
#' @param x An `hzr_stepwise` object.
#' @param ... Ignored.
#' @export
as.data.frame.hzr_stepwise <- function(x, ...) {
  x$steps
}


#' Extract the captured console trace from an `hzr_stepwise` fit
#'
#' Every run of [hzr_stepwise()] records the header, per-step lines,
#' and final summary regardless of the `trace` flag.  This accessor
#' returns the full character vector for display or logging.
#'
#' @param fit An `hzr_stepwise` object.
#' @return Character vector, one element per console line.
#' @export
stepwise_trace <- function(fit) {
  if (!inherits(fit, "hzr_stepwise")) {
    stop("`fit` must be an `hzr_stepwise` object.", call. = FALSE)
  }
  fit$trace_msg
}


# Local shape-count helper that tolerates multiphase as well as single
# distributions (the shared `.hzr_shape_parameter_count` in
# parity-helpers only handles the latter).
.hzr_stepwise_shape_count <- function(fit) {
  if (fit$spec$dist == "multiphase") {
    # For multiphase, "coefficient count" is sum of betas across phases
    # — i.e. theta length minus all non-beta slots.  Easiest path:
    # count columns in x_list.
    x_list <- fit$fit$x_list
    if (is.null(x_list)) return(length(fit$fit$theta))
    total_betas <- sum(vapply(x_list, function(m) {
      if (is.null(m)) 0L else ncol(m)
    }, integer(1L)))
    length(fit$fit$theta) - total_betas
  } else {
    .hzr_shape_parameter_count(fit$spec$dist,
                                control = fit$spec$control)
  }
}
