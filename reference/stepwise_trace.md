# Extract the captured console trace from an `hzr_stepwise` fit

Every run of
[`hzr_stepwise()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_stepwise.md)
records the header, per-step lines, and final summary regardless of the
`trace` flag. This accessor returns the full character vector for
display or logging.

## Usage

``` r
stepwise_trace(fit)
```

## Arguments

- fit:

  An `hzr_stepwise` object.

## Value

Character vector, one element per console line.
