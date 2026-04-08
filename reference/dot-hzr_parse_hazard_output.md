# Parse hazard binary output

Parses tabular output from the C hazard binary into a structured data
frame. Expects columns: parameter, estimate, StdErr (or SE), z-value (or
z), p-value (or pval). The parser is flexible and case-insensitive for
column matching.

## Usage

``` r
.hzr_parse_hazard_output(output_text, theta = NULL)
```

## Arguments

- output_text:

  Character vector (lines) from binary output file

- theta:

  Optional parameter vector for context/validation

## Value

Data frame with columns: parameter, estimate, std_err, z_stat, p_value
