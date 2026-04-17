# SAS Stepwise Parity Fixtures

This directory holds capture instructions and templates for
SAS-reference stepwise-selection runs.  The captured output lives as
`.rds` files under `inst/fixtures/` and is consumed by
`tests/testthat/test-stepwise-parity.R`.

Parity tests skip gracefully when a fixture is missing, so the package
still installs and passes CI without these files — they are only
required for claiming SAS-exact behaviour.

## What we capture

For each scenario (e.g. CABGKUL forward-Wald, AVC multiphase backward),
we want:

1. **The step-by-step trace** — which variable entered / dropped at
   each step, its Wald χ² statistic, degrees of freedom, p-value, and
   (for multiphase) which phase the move applied to.
2. **The final coefficient table** — estimates, standard errors,
   z-values, p-values at termination.
3. **Goodness-of-fit totals** — log-likelihood, AIC, iteration count.
4. **Meta** — SAS version, threshold settings, date captured, dataset,
   criterion, the exact `PROC HAZARD` invocation used.

Tolerances in the R-side comparison:

| Quantity                 | Tolerance           |
|--------------------------|---------------------|
| Sequence of (var, phase) | Exact               |
| Each step action         | Exact               |
| Wald χ² statistic        | 1e-3 (relative)     |
| p-value                  | 1e-3 (absolute)     |
| Final log-likelihood     | 1e-3 (absolute)     |
| Final coefficient        | 1e-2 (relative)     |

## Workflow

1. Run the SAS template in `cabgkul-forward-wald.sas` (adjust the
   `%LET` macros for other scenarios).
2. The template writes three text files into the SAS work library:
   - `stepwise_trace.csv` — one row per step
   - `stepwise_final.csv` — final coefficient table
   - `stepwise_meta.txt`  — key=value settings + GOF totals
3. In R, convert those three files into an `.rds` fixture with
   `TemporalHazard:::.hzr_build_stepwise_fixture()`, which
   validates the schema and writes to
   `inst/fixtures/stepwise-<scenario>.rds`.
4. Re-run `devtools::test()`; the parity test will pick up the
   fixture and assert each expectation.

See `schema.md` for the exact shape of the fixture list.
