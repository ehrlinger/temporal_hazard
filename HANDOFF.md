# Session Handoff — temporal_hazard

Last updated: 2026-04-10

## What Was Done This Session

### 1. hvtiPlotR: Kaplan-Meier Step Function Support

**Files changed:** `hvtiPlotR/R/hazard-plot.R`

Added `emp_geom` parameter to `hv_hazard()`, `plot.hv_hazard()`, and
`hazard_plot()`. When `emp_geom = "step"`, empirical KM data renders as
a proper step function (`geom_step()`) instead of discrete points
(`geom_point()`). Default is `"point"` for backward compat.

**Status:** Code complete, needs `devtools::install()` in hvtiPlotR
before temporal_hazard can use it.

### 2. temporal_hazard: Updated All KM Overlay Plots

Updated every KM-vs-parametric plot to use: - `emp_geom = "step"` for
proper KM step functions - Two distinct colors: blue (#0072B2) for
parametric, orange (#D55E00) for KM - Bottom-positioned legend with
`colour = NULL` (no legend title)

**Files changed:** - `R/hazard_api.R` — hazard() Rd example -
`vignettes/getting-started.qmd` -
`vignettes/hz-estimate-hazard-example.qmd` (AVC Weibull + KUL
Exponential) - `vignettes/hp-prediction-example.qmd`

**Status:** Code complete. To see the plots, reinstall hvtiPlotR first,
then reinstall temporal_hazard, then rebuild pkgdown.

### 3. Multiphase Implementation Plan

Created `MULTIPHASE_PLAN.md` with the full architecture for N-phase
hazard models using the generalized `decompos()` function from the
mixhazard repo.

**Key design decisions made:** - Prerelease — no backward compat
constraints - Every phase uses the same `decompos(t; t_half, nu, m)`
engine - N-phase (not hard-coded 3-phase) - Parameter names bridge SAS/C
conventions: `t_half`, `nu`, `m` - Additive cumulative hazard:
`H(t|x) = Σ μ_j(x) · Φ_j(t; θ_j)` - Phase types: “cdf” (early), “hazard”
(late), “constant” - Phase-specific covariates via formula argument in
[`hzr_phase()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase.md)

**Status:** Plan written, not yet implemented.

## What Needs Doing Next

### Immediate: Pre-Implementation Prep

1.  **Regenerate golden fixtures** — still pending from the Weibull
    optimizer reparameterization. Run:

    ``` r
    devtools::load_all()
    .hzr_create_synthetic_golden_fixtures()
    ```

    This will update the 3 .rds files in `inst/fixtures/` and fix the 4
    CI test failures.

2.  **Reinstall hvtiPlotR** — the step geom changes need to be
    installed:

    ``` r
    # From hvtiPlotR directory:
    devtools::install()
    ```

3.  **Rebuild pkgdown site** — to see the updated KM plots:

    ``` r
    # From temporal_hazard directory:
    devtools::install()
    pkgdown::build_site()
    ```

### Implementation: Multiphase (per MULTIPHASE_PLAN.md)

Build order (each step independently testable):

1.  **`R/decomposition.R`** — port `decompos()` from mixhazard as
    [`hzr_decompos()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_decompos.md).
    Add
    [`hzr_phase_cumhaz()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase_cumhaz.md)
    and
    [`hzr_phase_hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase_hazard.md).
    Test against mixhazard reference output.

2.  **`R/phase-spec.R`** —
    [`hzr_phase()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase.md)
    constructor. Simple S3 class.

3.  **`R/likelihood-multiphase.R`** — likelihood, gradient, multi-start
    optimizer. Follows the existing `likelihood-weibull.R` pattern.

4.  **`R/hazard_api.R`** — wire `dist = "multiphase"` dispatch with
    `phases` argument.

5.  **predict/summary** — add `decompose = TRUE` for per-phase output.

6.  **`R/argument_mapping.R`** — extend SAS/C → R mapping table.

7.  **Vignettes/examples** — decomposition plot showing phases adding
    together.

8.  **Parity tests** — validate against C binary output.

## Key Context

- The Weibull optimizer was reparameterized to use the C code’s (α, ψ,
  β) internally, with delta-method back-transform to user-facing (μ, ν,
  β). This is in `R/likelihood-weibull.R`.

- The `hazard` repo at `/Users/ehrlinj/Documents/GitHub/hazard/` has the
  original C source code for reference.

- The `mixhazard` repo at `/Users/ehrlinj/Documents/GitHub/mixhazard/`
  has the generalized `decompos()` function in `R/time_functions.R`.

- hvtiPlotR’s `hv_hazard()` handles all survival plotting. The package
  is at `/Users/ehrlinj/Documents/GitHub/hvtiPlotR/`.

## Files to Read First

When resuming:

    temporal_hazard/MULTIPHASE_PLAN.md     — full architecture plan
    temporal_hazard/R/likelihood-weibull.R — pattern for new likelihood files
    temporal_hazard/R/hazard_api.R         — main API to modify
    temporal_hazard/R/optimizer.R          — generic optimizer wrapper
    mixhazard/R/time_functions.R           — decompos() source to port
    hazard/src/utils/hzr_cum_haz_func.c   — original C phase functions
