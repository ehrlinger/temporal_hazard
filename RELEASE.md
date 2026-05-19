# Release process

Maintainer checklist for TemporalHazard. Not shipped (`.Rbuildignore`d).

## Versioning convention (the `.9000` strategy)

- **CRAN releases** carry a clean three-part version: `1.0.2`.
- **The development line** (everything between releases, on `main`) carries a
  fourth `.9000` component: `1.0.2.9000`. This makes it unambiguous whether a
  given install is the CRAN release or a dev build.
- A CRAN submission commit **always** strips `.9000` back to the clean release
  number. You never submit a `.9000` version to CRAN.

So the cycle is:

```
... 1.0.2.9000 (dev)  ->  1.0.3 (submit)  ->  1.0.3 accepted  ->  1.0.3.9000 (dev) ...
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

1. Tag the release and push the tag:

   ```sh
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```

2. Bump `main` to the next development version:

   ```r
   usethis::use_dev_version()   # DESCRIPTION -> X.Y.Z.9000
   ```

   or manually: set `DESCRIPTION` `Version:` to `X.Y.Z.9000` and add a NEWS
   heading:

   ```
   # TemporalHazard X.Y.Z.9000 (development version)
   ```

3. Commit to `main`: `chore: bump to X.Y.Z.9000 (dev)`.

All subsequent development happens at `X.Y.Z.9000` until the next release
strips it for submission.

## Historical note

Releases through `1.0.2` used clean three-part versions only; the `.9000`
dev convention is adopted starting **after `1.0.2` is accepted on CRAN**.
Existing tags: `v0.1.0`, `v0.9.3`, `v1.0.1` (so `v1.0.0` and `v1.0.2` should
be tagged at acceptance to keep the tag history consistent).
