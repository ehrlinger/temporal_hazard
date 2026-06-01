# Release process

Maintainer checklist for TemporalHazard. Not shipped (`.Rbuildignore`d).

## Versioning convention (the `.9000` strategy)

- **CRAN releases** carry a clean three-part version: `1.0.2`.
- **Development versions** carry a fourth `.9000` component (`1.0.2.9000`), so
  it is unambiguous whether a given install is the CRAN release or a dev build.
- A CRAN submission commit **always** strips `.9000` back to the clean release
  number. You never submit a `.9000` version to CRAN.

So the per-release cycle is:

```
... 1.0.2.9000 (dev)  ->  1.0.3 (submit)  ->  1.0.3 accepted  ->  1.0.3.9000 (dev) ...
```

See **Branch model** below for *which branch* carries which `.9000` line — this
package runs a two-branch split, so there are two: the released line on `main`
and the next-version line on `dev`.

## Branch model (dev / main split)

This is a CRAN package in an active, multi-PR development cycle, so it uses a
two-branch model. Not every package needs this — the split earns its keep only
when `main` has an external gate (CRAN) it must always satisfy. Internal /
no-CRAN packages should just work on `main` with feature branches and straight
semver. (Org-wide convention lives in the vault, `memory/r-package-branch-strategy.md`.)

| Branch | Version | Role |
|--------|---------|------|
| `main` | released `.9000` (e.g. `1.0.3.9000`) | the CRAN-released line; always shippable. Serves the live README, pkgdown site, and CRAN-facing docs. |
| `dev`  | next-version `.9000` (e.g. `1.1.0.9000`) | accumulates the next release; may hold breaking or half-finished work. |

**Why `.9000` on `main`?** `main` carries the released version plus any
between-release fixes (badges, README, figures). The `.9000` honestly marks
that `main` is *ahead of* the CRAN tarball, not the tarball itself. The pkgdown
navbar showing `1.0.3.9000` is expected and standard — the CRAN status badge
carries the released number.

**Routing:**

- **Shipped changes** (R code, `man/`, vignettes, tests) — branch off `dev`,
  PR to `dev`. They reach `main` only via a release.
- **Cosmetic / GitHub-facing fixes** (README prose, badges, figures, pkgdown
  config) — branch off `main`, PR to `main`. No version change beyond the
  existing `.9000`.
- **Release ops** (tag, GitHub Release) — directly on `main` (the documented
  exception to the branch+PR rule).

**Keep them in sync:** after a batch of `main`-side doc fixes, merge
`main -> dev` so they flow into the dev line. Skipping this is what produces
the "same fix, two parallel PRs" situation.

**Release flow:**

```
dev (1.1.0.9000) accumulates
  -> cut the suffix on dev: 1.1.0
  -> run the pre-submission gate (below), submit to CRAN
  -> on acceptance: merge dev -> main, tag v1.1.0, publish GitHub Release
  -> bump main -> 1.1.0.9000, open the next dev cycle
```

## Pre-submission checklist

1. `DESCRIPTION` `Version:` is a clean `X.Y.Z` (no `.9000`).
2. `NEWS.md` top heading is `# TemporalHazard X.Y.Z` and describes the
   user-visible changes since the last CRAN release.
3. `cran-comments.md` version heading matches `DESCRIPTION`; if a
   resubmission, each reviewer point is itemised with how it was addressed.
4. `devtools::document()` is clean and `man/` is in sync.
5. `R CMD check --as-cran` → 0 errors, 0 warnings, 0 notes
   (`devtools::check(args = "--as-cran")`).
6. `devtools::check_win_devel()` (and optionally `rhub::rhub_check()`).
   **This is the source of truth for the NOTE disposition** — the local
   `--as-cran` in step 5 does *not* run the CRAN incoming-feasibility /
   aspell step, so it will under-report. Reconcile the
   `## NOTE disposition` section of `cran-comments.md` against the
   *win-builder* `00check.log`, not the local result, before submitting.
   See "Known benign NOTE" below.
7. Submit: `devtools::submit_cran()` (writes `CRAN-SUBMISSION`).

### Known benign NOTE

Every submission produces one expected NOTE from the CRAN incoming check
that the local `--as-cran` does not show:

- **New submission** — until the package is accepted on CRAN. Drop this
  line from `cran-comments.md` once accepted.
- **Possibly misspelled words in DESCRIPTION** — proper nouns (`Naftel`,
  `Rajeswaran`), the acronym `UAB`, the `et al.` citation, and the domain
  term `multiphase`. All intentional; not misspellings. `inst/WORDLIST`
  feeds only the `spelling` test, **not** the CRAN aspell check, so this
  recurs by design — document it, don't try to suppress it.

The canonical wording lives in `cran-comments.md` `## NOTE disposition`;
keep the two in agreement.

## On CRAN acceptance

1. Tag the release and push the tag (annotated, mirroring existing tags):

   ```sh
   git tag -a vX.Y.Z -m "TemporalHazard vX.Y.Z — CRAN acceptance" <submitted-SHA>
   git push origin vX.Y.Z
   ```

   `<submitted-SHA>` is the `SHA:` line in `CRAN-SUBMISSION`.

2. Publish a GitHub Release from the tag, using that version's `NEWS.md`
   section as the body, marked as the latest release:

   ```sh
   gh release create vX.Y.Z --verify-tag --latest \
     --title "vX.Y.Z — <summary>" --notes-file <news-section.md>
   ```

3. **First CRAN release only — badge alignment.** Replace the manual
   shields `version` badge with the CRAN status badge and add the cranlogs
   download badges; the README `<!-- badges: start -->` block becomes:

   ```md
   <!-- badges: start -->
   [![CRAN status](https://www.r-pkg.org/badges/version/TemporalHazard)](https://CRAN.R-project.org/package=TemporalHazard)
   [![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/TemporalHazard)](https://CRAN.R-project.org/package=TemporalHazard)
   [![CRAN downloads (monthly)](https://cranlogs.r-pkg.org/badges/TemporalHazard)](https://CRAN.R-project.org/package=TemporalHazard)
   [![R-CMD-check](https://github.com/ehrlinger/temporal_hazard/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ehrlinger/temporal_hazard/actions/workflows/R-CMD-check.yaml)
   [![Codecov test coverage](https://codecov.io/gh/ehrlinger/temporal_hazard/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ehrlinger/temporal_hazard?branch=main)
   [![lint](https://github.com/ehrlinger/temporal_hazard/actions/workflows/lint.yaml/badge.svg)](https://github.com/ehrlinger/temporal_hazard/actions/workflows/lint.yaml)
   [![pkgdown site](https://img.shields.io/badge/docs-pkgdown-blue)](https://ehrlinger.github.io/temporal_hazard/)
   ![active](https://www.repostatus.org/badges/latest/active.svg)
   <!-- badges: end -->
   ```

   The manual version badge is now redundant (CRAN status badge is the
   single source of version truth), so also **delete the badge-sync CI**:

   ```sh
   git rm .github/workflows/update-version-badge.yaml
   ```

4. Bump `main` to the next development version:

   ```r
   usethis::use_dev_version()   # DESCRIPTION -> X.Y.Z.9000
   ```

   or manually: set `DESCRIPTION` `Version:` to `X.Y.Z.9000` and add a NEWS
   heading:

   ```
   # TemporalHazard X.Y.Z.9000 (development version)
   ```

5. Steps 3 and 4 change tracked files, so per the repo's mandatory git
   workflow they go through a branch + PR (not a direct `main` push):
   `chore: post-1.0.2 — CRAN badges + dev bump to X.Y.Z.9000`. Steps 1 and
   2 (tag push, GitHub Release) are the standard release-op exception.

All subsequent development happens at `X.Y.Z.9000` until the next release
strips it for submission.

## Historical note

Releases through `1.0.2` used clean three-part versions only; the `.9000`
dev convention is adopted starting **after `1.0.2` is accepted on CRAN**.
Existing tags: `v0.1.0`, `v0.9.3`, `v1.0.1` (so `v1.0.0` and `v1.0.2` should
be tagged at acceptance to keep the tag history consistent).
