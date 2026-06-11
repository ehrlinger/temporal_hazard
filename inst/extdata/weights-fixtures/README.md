# SAS Fractional-Weight Parity Fixture

This directory holds the capture payload for a SAS-reference run that
confirms `PROC HAZARD`'s `WEIGHT` statement matches R's
`hazard(weights = ...)` for **non-integer** (inverse-probability-style)
weights. It closes roadmap item **7a** / FIXTURE-GAP-LIST **B5**:

> Uniform non-integer weights (fractional sampling weights) — verify the
> SAS `WEIGHT` statement matches R behaviour for non-integer weights.

R-side correctness is already proven *without* SAS by the additive-split
and linear-scaling invariants in `tests/testthat/test-weights.R`. This
fixture is the external SAS confirmation. The parity test
(`tests/testthat/test-weights-sas-parity.R`) **skips gracefully when the
fixture is absent**, so the package installs and CI passes without any SAS
output — the fixture is only required to claim SAS-exact weighted behaviour.

## The model

A fixed (non-stepwise) two-phase fit on the AVC data:

- **Early** phase: Weibull CDF shape, `THALF`/`NU`/`M` fixed at the AVC
  unconditional values, covariate `AGE`.
- **Constant** phase: flat exponential hazard, covariate `AGE`.
- **Weight**: the fractional column `IPW` (`WEIGHT IPW;` in SAS,
  `weights = ipw` in R).
- Conservation of events on (`CONDITION=14`).

## The weight column

`payload/avc-weighted.csv` is the AVC analysis data with a deterministic,
non-integer `ipw` column appended. `ipw` is a closed-form function of
`age` only (logistic, centred at 60, scale 10 → range `(0.5, 1.5)`), with
**no RNG and no cohort-derived statistics**, so SAS and R read identical
values. Regenerate it with `make-avc-weighted-csv.R`. The exact formula is
illustrative; only the *shared, fractional* values matter for parity.

## Tolerances (R-side comparison)

| Quantity                       | Tolerance        |
|--------------------------------|------------------|
| Covariate coefficient (per phase) | 0.02 (relative) |
| Final log-likelihood           | 0.5 (absolute)   |

## Workflow

1. On a SAS box with the CCF HAZARD install, from `payload/`:
   ```sh
   ./run.sh           # auto-discovers HAZAPPS / HZ_MACROS; or set them
   ```
   This writes `payload/capture/avc-weighted.lst` and
   `payload/capture/weighted_meta.txt`.
2. Back in R, from the project root:
   ```r
   source("inst/extdata/weights-fixtures/parse-weighted-lst.R")
   ```
   This parses the listing's coefficient table + log-likelihood and calls
   `TemporalHazard:::.hzr_build_weighted_fixture()`, writing
   `inst/fixtures/weighted-avc-fractional.rds`.
3. Re-run `devtools::test(filter = "weights-sas-parity")`; the parity test
   picks up the fixture and asserts the covariate estimates and
   log-likelihood agree with SAS.

## Files

| File | Purpose |
|------|---------|
| `make-avc-weighted-csv.R` | Regenerates `payload/avc-weighted.csv` deterministically |
| `payload/avc-weighted.csv` | AVC analysis data + fractional `ipw` column |
| `payload/avc-weighted.sas` | `PROC HAZARD` weighted two-phase fit template |
| `payload/run.sh` | Discovers the HAZARD install and runs SAS |
| `payload/capture/` | Where the `.lst` + parsed outputs land |
| `parse-weighted-lst.R` | Parses the `.lst` and builds the `.rds` fixture |

The fixture list shape is defined by `.hzr_weighted_fixture_schema()` in
`R/weights-fixture.R`.
