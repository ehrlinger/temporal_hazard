# Finite-difference derivatives of G3 phase Phi and phi w.r.t. shape params

Computes the G3 cumulative intensity, its time derivative, and their
partial derivatives with respect to log_tau, gamma, alpha, and eta using
central finite differences.

## Usage

``` r
.hzr_g3_phase_derivatives(time, tau, gamma, alpha, eta, h = 1e-05)
```

## Arguments

- time:

  Numeric vector of positive times.

- tau:

  Positive scalar scale parameter.

- gamma:

  Positive scalar time exponent.

- alpha:

  Non-negative scalar shape parameter.

- eta:

  Positive scalar outer exponent.

- h:

  Relative step size for finite differences (default 1e-5).

## Value

Named list with Phi, phi, and 8 derivative vectors.
