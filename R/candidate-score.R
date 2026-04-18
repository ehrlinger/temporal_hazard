# candidate-score.R — Unified criterion wrapper for stepwise selection.
#
# Step 8.2 of STEPWISE-DESIGN.md.  Wraps the Wald helper (Phase 4b §8.1)
# and an AIC path behind a single signature so the forward and backward
# drivers can share structure.
#
# Score semantics
# ---------------
# Lower `score` is always "better" in the candidate's favour:
#
#   criterion="wald", mode="entry"  →  p_value(candidate for new coef)
#   criterion="wald", mode="drop"   →  1 − p_value(current)  (high p → easy drop)
#   criterion="aic",  mode="entry"  →  AIC(candidate) − AIC(current)
#   criterion="aic",  mode="drop"   →  ΔAIC_drop ≈ W − 2·df  (Wald→LR approx)
#
# For forward selection: argmin(score) picks the best candidate and the
# threshold test is `score < slentry` (Wald) or `score < 0` (AIC).
# For backward selection: argmin(score) picks the variable most worth
# dropping and the threshold test is `score < (1 − slstay)` (Wald) or
# `score < 0` (AIC).
#
# The list always carries both the raw `p_value` and `delta_aic` fields
# so the selection trace can show whichever is most informative
# regardless of the active criterion.

#' AIC for a fitted hazard model
#'
#' @param fit A fitted `hazard` object (`fit = TRUE`).
#' @return Numeric scalar AIC, or `NA_real_` if the log-likelihood is
#'   not available.
#'
#' @keywords internal
#' @noRd
.hzr_aic <- function(fit) {
  if (!inherits(fit, "hazard")) {
    stop("`fit` must be a `hazard` object.", call. = FALSE)
  }
  loglik <- fit$fit$objective
  if (is.null(loglik) || !is.finite(loglik)) {
    return(NA_real_)
  }

  # Count free parameters.  Multiphase fits with `fixed = "shapes"`
  # expose `fixed_mask` in fit_state, where TRUE marks a parameter held
  # fixed during optimisation (see `.hzr_optim_multiphase()`:
  # `best_result$fixed_mask <- !free_mask`).  Count only the FALSE
  # entries toward k; fall back to length(theta) when no mask is set.
  mask <- fit$fit$fixed_mask
  k <- if (!is.null(mask) && length(mask) == length(fit$fit$theta)) {
    sum(!as.logical(mask))  # fixed_mask element TRUE == fixed
  } else {
    length(fit$fit$theta)
  }

  -2 * loglik + 2 * k
}


#' Score a stepwise candidate under Wald or AIC criterion
#'
#' Produces a unified "smaller is better" score for either adding new
#' coefficients to a model (`mode = "entry"`) or dropping coefficients
#' from the current model (`mode = "drop"`), under a Wald or AIC
#' criterion.
#'
#' @param criterion Either `"wald"` or `"aic"`.
#' @param mode Either `"entry"` or `"drop"`.
#' @param current Fitted `hazard` object that is the starting point for
#'   the step.  Always required (AIC baseline; Wald target for drops).
#' @param candidate Fitted `hazard` object containing the proposed new
#'   coefficient(s).  Required when `mode = "entry"`; ignored otherwise.
#' @param names Character vector of coefficient names being tested.
#'   Must exist in the target fit (`candidate` for entries, `current`
#'   for drops).
#'
#' @return A list with:
#' \describe{
#'   \item{score}{Numeric, smaller is better (see file header).}
#'   \item{criterion}{Echoed input.}
#'   \item{mode}{Echoed input.}
#'   \item{stat}{Wald z (scalar) or χ² (joint).}
#'   \item{df}{Integer degrees of freedom.}
#'   \item{p_value}{Always populated when computable.}
#'   \item{delta_aic}{Always populated when computable.}
#'   \item{names}{Echoed input.}
#' }
#'
#' @keywords internal
#' @noRd
.hzr_candidate_score <- function(criterion = c("wald", "aic"),
                                  mode      = c("entry", "drop"),
                                  current,
                                  candidate = NULL,
                                  names) {
  criterion <- match.arg(criterion)
  mode      <- match.arg(mode)

  if (!inherits(current, "hazard")) {
    stop("`current` must be a fitted `hazard` object.", call. = FALSE)
  }
  if (mode == "entry") {
    if (is.null(candidate) || !inherits(candidate, "hazard")) {
      stop("`candidate` must be a fitted `hazard` object when mode = 'entry'.",
           call. = FALSE)
    }
    test_fit <- candidate
  } else {
    test_fit <- current
  }

  wald <- .hzr_wald_p(test_fit, names)

  # Wald → LR approximation:  W = z² when df = 1;
  # otherwise W is already the χ² statistic.
  W <- if (wald$df == 1L) {
    if (is.finite(wald$stat)) wald$stat^2 else NA_real_
  } else {
    if (is.finite(wald$stat)) wald$stat else NA_real_
  }

  # AIC components
  if (criterion == "aic" && mode == "entry") {
    aic_cur <- .hzr_aic(current)
    aic_can <- .hzr_aic(candidate)
    delta   <- if (is.finite(aic_cur) && is.finite(aic_can)) {
      aic_can - aic_cur
    } else {
      NA_real_
    }
  } else if (criterion == "aic" && mode == "drop") {
    delta <- if (is.finite(W)) W - 2 * wald$df else NA_real_
  } else if (mode == "drop") {
    # Wald + drop: compute ΔAIC approx anyway so the trace can show it.
    delta <- if (is.finite(W)) W - 2 * wald$df else NA_real_
  } else {
    # Wald + entry: no ΔAIC available without refitting;  leave as NA.
    delta <- NA_real_
  }

  score <- if (criterion == "wald") {
    if (mode == "entry") {
      wald$p_value           # smaller p = stronger entry case
    } else {
      if (is.na(wald$p_value)) NA_real_ else 1 - wald$p_value
    }                        # larger p = easier drop, so 1-p is small
  } else {
    delta                    # AIC: smaller ΔAIC = better outcome
  }

  list(
    score     = score,
    criterion = criterion,
    mode      = mode,
    stat      = wald$stat,
    df        = wald$df,
    p_value   = wald$p_value,
    delta_aic = delta,
    names     = names
  )
}
