# Clamp probabilities away from 0 and 1

Clamp probabilities away from 0 and 1

## Usage

``` r
hzr_clamp_prob(p, eps = 1e-12)
```

## Arguments

- p:

  Numeric vector of probabilities.

- eps:

  Small positive tolerance.

## Value

Numeric vector bounded to `[eps, 1 - eps]`.
