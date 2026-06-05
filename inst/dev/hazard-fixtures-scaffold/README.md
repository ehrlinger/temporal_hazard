# hazard-fixtures

Private test fixture repository for cross-codebase parity between the
[hazard](https://github.com/ehrlinger/hazard) C binary and the
[TemporalHazard](https://github.com/ehrlinger/temporal_hazard) R package.

---

## PHI policy — read before adding anything

**Allowed in this repo:**
- `.lst` files — SAS HAZARD output captures (aggregate statistics: log-likelihood,
  parameter estimates, standard errors, variance-covariance matrix, observation counts)
- `.sas` files — SAS job scripts (code only, no inline data)
- `.rds` files — R fixture objects derived from the above
- This README and supporting documentation

**Never commit:**
- Raw patient-level flat files (e.g., `omc`, `avc`, `kul`)
- SAS datasets (`.sas7bdat`, `.ssd01`)
- Any file containing one row per patient
- Any file traceable to an individual

The R datasets shipped with TemporalHazard (`cabgkul`, `avc`, `valves`, `tga`, `omc`)
are de-identified and public (on CRAN). They are used directly from the installed
package — they do not live here.

---

## Repository layout

```
hazard-fixtures/
├── README.md
├── .gitignore              # enforces PHI policy at the file level
├── examples/               # .lst + .sas pairs, one per SAS job
│   ├── hz.deadp.KUL.sas
│   ├── hz.deadp.KUL.lst
│   ├── hz.death.AVC.sas
│   ├── hz.death.AVC.lst
│   ├── hz.te123.OMC.sas    # .lst IS committed (aggregate); re-running needs raw omc data
│   ├── hz.te123.OMC.lst
│   ├── hz.tm123.OMC.sas
│   ├── hz.tm123.OMC.lst
│   ├── hm.deadp.VALVES.sas
│   ├── hm.deadp.VALVES.lst
│   ├── hm.death.patient.sas
│   ├── hm.death.patient.lst
│   ├── hm.death.AVC.deciles.sas
│   ├── hm.death.AVC.deciles.lst
│   ├── hm.death.AVC.sas
│   ├── hm.death.AVC.lst
│   ├── ac.death.AVC.sas
│   ├── ac.death.AVC.lst
│   ├── hp.death.AVC.sas
│   ├── hp.death.AVC.lst
│   ├── hp.death.AVC.hm1.sas
│   ├── hp.death.AVC.hm1.lst
│   ├── hp.death.AVC.hm2.sas
│   ├── hp.death.AVC.hm2.lst
│   ├── hs.death.AVC.hm1.sas
│   ├── hs.death.AVC.hm1.lst
│   ├── bs.death.AVC.sas
│   └── bs.death.AVC.lst
├── r-fixtures/             # .rds golden files for R-only parity tests
│   └── stepwise-cabgkul-forward-wald.rds
└── .github/
    └── workflows/
        └── validate-fixtures.yaml
```

---

## Environment variables

| Variable | Points to | Used by |
|---|---|---|
| `HAZARD_EXAMPLES_DIR` | `<repo>/examples/` | R `test-sas-parity.R`, C `validate_corpus.sh` |
| `HAZARD_OMC_RAW` | raw `omc` flat file on CCF network share | R `test-sas-parity.R` (PRIMISOL tests only) |
| `TEMPORAL_HAZARD_BIN` | compiled `hazard` binary | R `test-parity-c-binary.R` |

`HAZARD_OMC_RAW` is **never set in CI** — the `hz.te123.OMC` and `hz.tm123.OMC`
parity tests skip automatically when it is unset.

### Setting up locally (macOS/Linux)

Add to your shell profile or `.Renviron`:

```bash
export HAZARD_EXAMPLES_DIR="$HOME/Documents/GitHub/hazard-fixtures/examples"
export HAZARD_OMC_RAW="/path/to/ccf/network/share/omc"   # local machine only — never set in CI or commit the path
export TEMPORAL_HAZARD_BIN="$HOME/Documents/GitHub/hazard/hazard"
```

---

## Capture vintage and refresh policy

Each `.lst` file was captured from a specific binary version on a specific platform.
The capture vintage is recorded in the file header comment and in the table below.

| Fixture group | Binary version | Platform | Captured |
|---|---|---|---|
| All primary fits | v4.4.6 | macOS arm64 | 2026-05-12 |
| Windows reference | v4.4.5 | Windows x64 + SAS | 2026-04-28 |

When a new binary version changes output format or numerical results, re-capture
on the same platform and update this table. Note: the validate-fixtures CI only
checks that each `.lst` still **parses** — it does not compare numeric results
across versions, so it will catch format-breaking changes but not silent
numerical drift. Cross-version numerical regression is the consuming repos'
parity tests' job (run with the env vars above).

---

## Adding a new fixture

1. Write or locate the `.sas` job script. Confirm it contains no inline data.
2. Run it against the C binary locally with access to the appropriate dataset.
3. Copy the `.lst` output into `examples/`. Verify it parses cleanly:
   ```r
   source("tests/testthat/helper-sas-parity.R")  # from temporal_hazard repo
   .hzr_parse_sas_lst("examples/my.new.job.lst")
   ```
4. Commit `.sas` + `.lst`. Do not commit any data files.
5. Write the corresponding `test_that()` block in `temporal_hazard/tests/testthat/test-sas-parity.R`.

---

## CI access for public repos

Both `hazard` and `temporal_hazard` have a `.github/workflows/fixture-parity.yaml`
that checks out this repo using a deploy key. See the setup section in each repo's
workflow file for key generation instructions. The deploy key grants read-only access
and is stored as the `FIXTURE_REPO_DEPLOY_KEY` secret on each public repo.
