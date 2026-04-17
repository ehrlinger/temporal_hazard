# Stepwise Fixture Schema

Each fixture `.rds` file deserialises to a named list with four
components: `meta`, `steps`, `final`, and `scope`.  All fields are
required unless marked optional.

```r
list(
  meta = list(
    dataset       = "cabgkul",          # name of the data frame used
    dist          = "weibull",          # or "exponential", "multiphase", ...
    criterion     = "wald",             # or "aic"
    direction     = "forward",          # or "backward", "both"
    slentry       = 0.30,               # entry threshold (criterion = wald)
    slstay        = 0.20,               # retention threshold (criterion = wald)
    sas_version   = "9.4",
    captured_on   = "2026-04-17",       # ISO-8601
    proc_hazard   = "PROC HAZARD ...",  # exact SAS invocation for audit
    notes         = NULL                # optional free text
  ),

  scope = list(
    # For single-distribution models: character vector of candidate
    # names.  For multiphase: named list keyed by phase.
    candidates = c("age", "mal", "shock", "nyha", "inc_surg"),
    force_in   = character(),
    force_out  = character()
  ),

  steps = data.frame(
    step_num  = integer(),     # 1..N
    action    = character(),   # "enter" | "drop"
    variable  = character(),
    phase     = character(),   # NA_character_ for single-dist
    stat      = numeric(),     # Wald chi-square
    df        = integer(),     # degrees of freedom (usually 1)
    p_value   = numeric(),     # two-sided
    stringsAsFactors = FALSE
  ),

  final = list(
    coef = data.frame(
      variable  = character(),  # same order as final model's print
      phase     = character(),  # NA_character_ for single-dist
      estimate  = numeric(),
      std_error = numeric(),
      z_stat    = numeric(),
      p_value   = numeric(),
      stringsAsFactors = FALSE
    ),
    logLik     = numeric(1L),   # at termination
    aic        = numeric(1L),   # optional; computed if absent
    iterations = integer(1L)    # SAS's iteration count (reported)
  )
)
```

## Notes

- SAS `PROC HAZARD` reports Wald chi-square; we store it directly
  rather than converting to z.  The R helper squares z for df = 1 to
  compare.
- For multiphase scenarios, `phase` must contain the SAS phase name
  as it appears in the output (e.g. `"EARLY"`, `"LATE"`, `"CONSTANT"`)
  — the R-side test lowercases and matches to `names(fit$spec$phases)`.
- `final$coef` must be in the same order as the R fit reports it when
  all entered variables are present.  The comparator sorts by name
  before checking, so ordering is for readability only.
- Tolerances live in `test-stepwise-parity.R`; see
  `README.md` for the current values.
