# Release process

Maintainer checklist for TemporalHazard. Not shipped (`.Rbuildignore`d).

## Versioning convention (the `.9000` strategy)

- **CRAN releases** carry a clean three-part version: `1.0.2`.
- **The development line** (everything between releases, on `main`)
  carries a fourth `.9000` component: `1.0.2.9000`. This makes it
  unambiguous whether a given install is the CRAN release or a dev
  build.
- A CRAN submission commit **always** strips `.9000` back to the clean
  release number. You never submit a `.9000` version to CRAN.

So the cycle is:

    ... 1.0.2.9000 (dev)  ->  1.0.3 (submit)  ->  1.0.3 accepted  ->  1.0.3.9000 (dev) ...

## Pre-submission checklist

1.  `DESCRIPTION` `Version:` is a clean `X.Y.Z` (no `.9000`).
2.  `NEWS.md` top heading is `# TemporalHazard X.Y.Z` and describes the
    user-visible changes since the last CRAN release.
3.  `cran-comments.md` version heading matches `DESCRIPTION`; if a
    resubmission, each reviewer point is itemised with how it was
    addressed; test environments and NOTE disposition reflect reality.
4.  `devtools::document()` is clean and `man/` is in sync.
5.  `R CMD check --as-cran` → 0 errors, 0 warnings, 0 notes
    (`devtools::check(args = "--as-cran")`).
6.  Optional wider checks: `devtools::check_win_devel()`,
    `rhub::rhub_check()`.
7.  Submit: `devtools::submit_cran()` (writes `CRAN-SUBMISSION`).

## On CRAN acceptance

1.  Tag the release and push the tag (annotated, mirroring existing
    tags):

    ``` sh
    git tag -a vX.Y.Z -m "TemporalHazard vX.Y.Z — CRAN acceptance" <submitted-SHA>
    git push origin vX.Y.Z
    ```

    `<submitted-SHA>` is the `SHA:` line in `CRAN-SUBMISSION`.

2.  Publish a GitHub Release from the tag, using that version’s
    `NEWS.md` section as the body, marked as the latest release:

    ``` sh
    gh release create vX.Y.Z --verify-tag --latest \
      --title "vX.Y.Z — <summary>" --notes-file <news-section.md>
    ```

3.  **First CRAN release only — badge alignment.** Replace the manual
    shields `version` badge with the CRAN status badge and add the
    cranlogs download badges; the README `<!-- badges: start -->` block
    becomes:

    ``` md
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
    single source of version truth), so also **delete the badge-sync
    CI**:

    ``` sh
    git rm .github/workflows/update-version-badge.yaml
    ```

4.  Bump `main` to the next development version:

    ``` r

    usethis::use_dev_version()   # DESCRIPTION -> X.Y.Z.9000
    ```

    or manually: set `DESCRIPTION` `Version:` to `X.Y.Z.9000` and add a
    NEWS heading:

        # TemporalHazard X.Y.Z.9000 (development version)

5.  Steps 3 and 4 change tracked files, so per the repo’s mandatory git
    workflow they go through a branch + PR (not a direct `main` push):
    `chore: post-1.0.2 — CRAN badges + dev bump to X.Y.Z.9000`. Steps 1
    and 2 (tag push, GitHub Release) are the standard release-op
    exception.

All subsequent development happens at `X.Y.Z.9000` until the next
release strips it for submission.

## Historical note

Releases through `1.0.2` used clean three-part versions only; the
`.9000` dev convention is adopted starting **after `1.0.2` is accepted
on CRAN**. Existing tags: `v0.1.0`, `v0.9.3`, `v1.0.1` (so `v1.0.0` and
`v1.0.2` should be tagged at acceptance to keep the tag history
consistent).
