# Branch Protection Guidance

This repository should enforce branch protection on `main` (and `master` if still used) so pull requests cannot merge unless CI passes.

## Recommended Rules

Apply these settings in GitHub:

1. Settings -> Branches -> Add branch protection rule.
2. Branch name pattern: `main`.
3. Enable **Require a pull request before merging**.
4. Enable **Require approvals** (recommended: 1 minimum).
5. Enable **Dismiss stale pull request approvals when new commits are pushed**.
6. Enable **Require status checks to pass before merging**.
7. Enable **Require branches to be up to date before merging**.
8. Enable **Require conversation resolution before merging**.
9. Optional: Enable **Include administrators** for strict enforcement.

## Required Status Checks

Mark the following checks as required for PRs into `main`:

- `lint`
- `test-coverage`
- `build-and-deploy` (pkgdown workflow PR validation)
- `check-badge-sync` (README version badge consistency)
- All `R-CMD-check` matrix checks (each appears as a separate check run)

If check names vary by matrix label, require each visible `R-CMD-check (...)` entry.

## Why this matters

This ensures failures are caught on PR push/update rather than after merge to `main`, including pkgdown index drift and documentation consistency issues.
