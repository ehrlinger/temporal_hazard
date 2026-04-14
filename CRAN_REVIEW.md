# CRAN Readiness Review (2026-04-13)

Repository: temporal_hazard  
Branch reviewed: feature/analytic-gradients  
Reviewer: GitHub Copilot (GPT-5.3-Codex)

## Scope

Thorough review for CRAN readiness with emphasis on:
- consistency and R coding practices
- test coverage reasonableness
- documentation completeness and consistency
- CRAN check behavior from a built source tarball

## Validation commands run

- `R CMD build .`
- `R CMD check --as-cran --no-manual TemporalHazard_0.9.1.tar.gz`
- `Rscript -e 'testthat::test_local(".", reporter = "summary")'`
- `Rscript -e 'x <- covr::package_coverage(); covr::percent_coverage(x)'`
- `Rscript -e 'lintr::lint_package()'`

## Executive summary

- Built-tarball check is close to CRAN-ready: **0 ERROR, 0 WARNING, 1 NOTE**.
- The remaining NOTE is CRAN incoming/new-submission context plus a redirected URL in `inst/CITATION`.
- Main code-quality issue found: documentation default for `fit` is inconsistent with implementation.
- Test breadth is strong, but measured line coverage is moderate for a numerics-heavy package and can be improved in specific areas.

## Findings (prioritized)

### High

1. Public documentation contradicts function default for `fit`.
- Implementation default is `fit = FALSE` in `R/hazard_api.R`.
- Docs say default TRUE in both roxygen and Rd.
- Impact: user confusion and behavior mismatch in examples and downstream usage.
- Evidence:
  - `R/hazard_api.R` (function signature uses `fit = FALSE`)
  - `R/hazard_api.R` (`@param fit` text says default TRUE)
  - `man/hazard.Rd` (`fit` item says default TRUE)

### Medium

2. Multiphase optimizer uses non-deterministic reseeding inside optimization loop.
- `set.seed(NULL)` appears in multi-start path.
- Impact: run-to-run variability can hinder reproducibility for statistical workflows and debugging.
- Evidence:
  - `R/likelihood-multiphase.R` (multi-start loop with `set.seed(NULL)`)

3. Coverage level is moderate for a statistical core package.
- Observed total coverage around **62.4% to 62.8%**.
- Lowest coverage includes helper/core internals (for example, parts of multiphase internals and formula helpers), and fixture-generation code.
- Impact: increased regression risk in edge-case numerical behavior.
- Evidence from covr runs:
  - low file coverage reported for `R/golden_fixtures.R`, `R/formula-helpers.R`, `R/likelihood-multiphase.R`, `R/decomposition.R`

### Low

4. Vignette metadata style is inconsistent across `.qmd` files.
- `vignettes/ar-architecture.qmd` uses YAML `vignette:` metadata.
- Other vignettes rely on commented `%\Vignette...` blocks instead of YAML `vignette:`.
- Impact: maintenance and build-convention inconsistency.
- Evidence:
  - `vignettes/ar-architecture.qmd`
  - `vignettes/getting-started.qmd`
  - `vignettes/fitting-hazard-models.qmd`
  - `vignettes/prediction-visualization.qmd`
  - `vignettes/inference-diagnostics.qmd`
  - `vignettes/mf-mathematical-foundations.qmd`
  - `vignettes/sas-to-r-migration.qmd`

5. Remaining CRAN NOTE includes redirected URL in citation metadata.
- URL in citation redirects (HTTP 200 after move).
- Impact: minor incoming NOTE; usually easy to clean up.
- Evidence:
  - `inst/CITATION`
  - `TemporalHazard.Rcheck/00check.log`

## What looks good

- Built tarball check is clean apart from one NOTE.
- Namespace/documentation checks passed.
- Test suite breadth is strong: multiple distributions, censoring types, multiphase behavior, parity-focused tests.
- Lintr run found no lint issues.

## Recommended follow-up tasks

1. Fix docs/behavior mismatch for `fit` default.
- Choose one source of truth:
  - either change docs to default FALSE
  - or change API default to TRUE (if intended)
- Regenerate docs (`roxygen2`) and re-check.

2. Improve optimizer reproducibility policy.
- Consider deterministic multi-start when a seed is set by user.
- Optionally expose a control parameter for RNG seed policy.

3. Raise coverage on high-value logic.
- Prioritize exported method behavior (`vcov.hazard`, `coef.hazard`) and critical multiphase edge paths.
- Add targeted tests for numerical edge cases in currently lower-coverage files.

4. Normalize vignette metadata approach.
- Use consistent Quarto vignette YAML across all vignette `.qmd` files.

5. Remove/resolve remaining incoming NOTE where practical.
- Update citation URL to canonical non-redirecting target.

## Suggested release gate before CRAN submit

- `R CMD build .`
- `R CMD check --as-cran TemporalHazard_<version>.tar.gz`
- Ensure result is 0 ERROR / 0 WARNING and only expected NOTE(s) for new submission.
- Re-run tests and (optionally) coverage report after fixes.
