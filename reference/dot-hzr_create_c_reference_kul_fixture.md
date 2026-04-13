# Create golden fixture for C binary parity (KUL CABG dataset)

Loads the KUL cardiac surgery dataset from inst/extdata/cabgkul.csv and
stores it alongside the C HAZARD binary reference output for parity
testing. This fixture does NOT fit a model — it stores the data and C
reference values so that tests can evaluate the R log-likelihood at
C-estimated parameters and/or fit the model and compare.

## Usage

``` r
.hzr_create_c_reference_kul_fixture(output_dir = NULL)
```
