# Migrating from SAS HAZARD to TemporalHazard

``` r
library(TemporalHazard)
```

## Overview

The legacy HAZARD program is a SAS macro (`%HAZARD(...)`) that wraps a
compiled C executable implementing the multiphase parametric hazard
model of Blackstone, Naftel, and Turner (1986). This vignette is the
canonical reference for translating a SAS HAZARD analysis into native R
using `TemporalHazard`.

The full formal argument mapping table is available programmatically:

``` r
knitr::kable(
  TemporalHazard::hzr_argument_mapping(),
  caption = "Formal argument map: SAS HAZARD/C → hazard()",
  col.names = c(
    "SAS Statement", "Legacy Input", "C Concept",
    "R Parameter", "Required", "Expected Type",
    "Transform Rule", "Status", "Notes"
  )
)
```

| SAS Statement | Legacy Input               | C Concept                   | R Parameter                   | Required | Expected Type                | Transform Rule                                                                         | Status      | Notes                                                                                                    |
|:--------------|:---------------------------|:----------------------------|:------------------------------|:---------|:-----------------------------|:---------------------------------------------------------------------------------------|:------------|:---------------------------------------------------------------------------------------------------------|
| HAZARD        | TIME variable              | obs time array              | time                          | TRUE     | numeric vector               | pass through as numeric                                                                | implemented | Core observation time input.                                                                             |
| HAZARD        | EVENT/censor variable      | event indicator array       | status                        | TRUE     | numeric/logical vector       | coerce to numeric 0/1                                                                  | implemented | Event indicator currently retained as numeric in object$data$status.                                     |
| HAZARD        | X covariate block          | design matrix               | x                             | FALSE    | numeric matrix or data.frame | data.frame -\> data.matrix                                                             | implemented | Future versions will support richer design encoding helpers.                                             |
| HAZARD        | initial parameters         | parameter vector            | theta                         | FALSE    | numeric vector               | length must equal ncol(x) when x is present                                            | implemented | Used by predict.hazard as coefficient vector.                                                            |
| HAZARD        | baseline distribution      | phase distribution selector | dist                          | FALSE    | character scalar             | normalized lower-case label                                                            | implemented | Current default is ‘weibull’; more options planned.                                                      |
| HAZARD        | control options            | optimizer/control struct    | control                       | FALSE    | named list                   | stored in spec\$control                                                                | implemented | Control list is stored and reserved for optimizer parity.                                                |
| HAZARD        | additional legacy options  | misc legacy switches        | …                             | FALSE    | named arguments              | stored in legacy_args for parity                                                       | implemented | Supports legacy-style pass-through options during migration.                                             |
| TIME          | t                          | time vector                 | time                          | TRUE     | numeric vector               | pass through                                                                           | implemented | Canonical SAS migration uses TIME= mapping.                                                              |
| EVENT         | status                     | event vector                | status                        | TRUE     | numeric/logical vector       | coerce to numeric                                                                      | implemented | Canonical SAS migration uses EVENT= mapping.                                                             |
| PARMS         | theta0                     | starting coef               | theta                         | FALSE    | numeric vector               | map PARMS/INITIAL to theta                                                             | planned     | SAS PARMS syntax parser not yet implemented.                                                             |
| DIST          | dist                       | dist selector               | dist                          | FALSE    | character scalar             | map DIST= to dist                                                                      | implemented | SAS DIST keyword maps directly to dist.                                                                  |
| HAZARD        | phases (3-phase structure) | 3-phase Early/Const/Late    | phases (list of hzr_phase())  | FALSE    | list of hzr_phase objects    | list(early=hzr_phase(‘cdf’,…), constant=hzr_phase(‘constant’), late=hzr_phase(‘g3’,…)) | implemented | Use dist=‘multiphase’ with phases argument. N-phase generalization of legacy 3-phase model.              |
| HAZARD        | MU_1, MU_2, MU_3           | per-phase scale factors     | mu (via exp(log_mu) in theta) | FALSE    | numeric (per-phase)          | exp(alpha_j) in internal parameterization; estimated on log scale                      | implemented | Each phase has its own scale mu_j(x) = exp(alpha_j + x\*beta_j). Starting value via hzr_phase().         |
| G1            | THALF / RHO (early)        | early half-life             | hzr_phase(t_half=)            | FALSE    | positive scalar              | maps directly to hzr_phase(t_half=) starting value                                     | implemented | Half-life: time at which G(t_half) = 0.5. Same concept as SAS RHO/THALF.                                 |
| G1            | NU (early)                 | early time exponent         | hzr_phase(nu=)                | FALSE    | numeric scalar               | maps directly to hzr_phase(nu=) starting value                                         | implemented | Time exponent controlling rate dynamics. Same parameter name as SAS early NU.                            |
| G1            | M (early)                  | early shape                 | hzr_phase(m=)                 | FALSE    | numeric scalar               | maps directly to hzr_phase(m=) starting value                                          | implemented | Shape exponent controlling distributional form. Same parameter name as SAS early M.                      |
| G1            | DELTA (early)              | early time transform        | (absorbed by decompos)        | FALSE    | numeric scalar               | time transform B(t) = (exp(delta\*t)-1)/delta absorbed into decompos shape             | implemented | The C DELTA controlled B(t) = (exp(delta\*t)-1)/delta. This transform is absorbed by decompos().         |
| G2            | G2 constant phase          | constant hazard rate phase  | hzr_phase(‘constant’)         | FALSE    | hzr_phase(‘constant’)        | hzr_phase(‘constant’) with no shape parameters                                         | implemented | Flat background rate. No shape parameters estimated. SAS G2 equivalent.                                  |
| G3            | TAU (late)                 | late G3 scale               | hzr_phase(‘g3’, tau=)         | FALSE    | positive scalar              | maps directly to hzr_phase(‘g3’, tau=) for late phase                                  | implemented | Late-phase G3 scale parameter. Maps directly to hzr_phase(‘g3’, tau=).                                   |
| G3            | GAMMA (late)               | late G3 time exponent       | hzr_phase(‘g3’, gamma=)       | FALSE    | numeric scalar               | maps directly to hzr_phase(‘g3’, gamma=) for late phase                                | implemented | Late-phase G3 time exponent. Maps directly to hzr_phase(‘g3’, gamma=).                                   |
| G3            | ALPHA (late)               | late G3 shape               | hzr_phase(‘g3’, alpha=)       | FALSE    | numeric scalar               | maps directly to hzr_phase(‘g3’, alpha=) for late phase                                | implemented | Late-phase G3 shape parameter. alpha=0 gives exponential case. Maps directly to hzr_phase(‘g3’, alpha=). |
| G3            | ETA (late)                 | late G3 outer exponent      | hzr_phase(‘g3’, eta=)         | FALSE    | numeric scalar               | maps directly to hzr_phase(‘g3’, eta=) for late phase                                  | implemented | Late-phase G3 outer exponent. Maps directly to hzr_phase(‘g3’, eta=).                                    |

Formal argument map: SAS HAZARD/C → hazard()

------------------------------------------------------------------------

## Statement-by-statement mapping

### `PROC HAZARD DATA=`

``` sas
PROC HAZARD DATA=AVCS NOCOV NOCOR CONDITION=14;
```

Maps to
[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
`control` list:

``` r
fit <- hazard(
  ...,
  control = list(
    nocov      = TRUE,   # suppress covariance output
    nocor      = TRUE,   # suppress correlation output
    condition  = 14      # CONDITION= switch
  )
)
```

Additional `PROC HAZARD` options with no direct R equivalent yet are
passed through `...` as named arguments and stored in `fit$legacy_args`.

------------------------------------------------------------------------

### `TIME`

``` sas
TIME INT_DEAD;
```

Maps directly to `time`:

``` r
fit <- hazard(
  time   = avcs$INT_DEAD,
  ...
)
```

- Must be a non-negative numeric vector in the same time unit as the
  original analysis (months, in the AVC example).
- Missing values are not permitted.

------------------------------------------------------------------------

### `EVENT`

``` sas
EVENT DEAD;
```

Maps to `status`:

``` r
fit <- hazard(
  status = avcs$DEAD,   # 1 = event, 0 = censored
  ...
)
```

- `TemporalHazard` coerces `status` to `numeric`; logical vectors are
  also accepted.
- The HAZARD convention uses `1` = occurred, `0` = censored/alive at
  last contact.

------------------------------------------------------------------------

### `PARMS`

``` sas
PARMS MUE=0.3504743 THALF=0.1905077 NU=1.437416 M=1 FIXM
      MUC=4.391673E-07;
```

Maps to `theta` (coefficient/parameter vector) and `control` (fix
flags):

``` r
fit <- hazard(
  theta   = c(MUE = 0.3504743, THALF = 0.1905077, NU = 1.437416,
              M = 1,           MUC   = 4.391673e-07),
  control = list(
    fix = c("M")   # FIXM → freeze M during optimization
  ),
  ...
)
```

> **Note:** Full SAS `PARMS` syntax is not mirrored one-for-one in the
> public API. In the current package, supply the parameter vector
> directly via `theta` and record fixed-parameter intent in
> `control$fix`. The legacy parity helpers do generate SAS-style control
> text when comparing against the historical binaries.

------------------------------------------------------------------------

### `EARLY` and `CONSTANT` covariate blocks

``` sas
EARLY  AGE=-0.03205774, COM_IV=1.336675, MAL=0.6872028,
       OPMOS=-0.01963377, OP_AGE=0.0002086689, STATUS=0.5169533;
CONSTANT INC_SURG=1.375285, ORIFICE=3.11765, STATUS=1.054988;
```

In HAZARD, `EARLY` and `CONSTANT` define phase-specific covariate
coefficients. In `TemporalHazard` these are unified into a single design
matrix `x` and coefficient vector `theta` during M1. Phase assignment
will be formalised in M2 when the multi-phase likelihood is implemented.

**Current convention (M1):** combine all covariates into `x` and supply
the corresponding starting coefficients in `theta`.

``` r
# Build design matrix from the AVC data set
X <- data.matrix(avcs[, c("AGE", "COM_IV", "MAL", "OPMOS", "OP_AGE",
                           "STATUS", "INC_SURG", "ORIFICE")])

# Starting values from SAS EARLY + CONSTANT blocks combined
theta_start <- c(
  AGE      = -0.03205774,
  COM_IV   =  1.336675,
  MAL      =  0.6872028,
  OPMOS    = -0.01963377,
  OP_AGE   =  0.0002086689,
  STATUS   =  0.5169533,   # EARLY phase coefficient
  INC_SURG =  1.375285,
  ORIFICE  =  3.11765
)

fit <- hazard(
  time   = avcs$INT_DEAD,
  status = avcs$DEAD,
  x      = X,
  theta  = theta_start,
  dist   = "weibull"
)
```

------------------------------------------------------------------------

### `SELECTION`

``` sas
SELECTION SLE=0.2 SLS=0.1;
```

Stepwise variable selection is **not implemented** in `TemporalHazard`.
For now, resolve the final variable set outside the package and then fit
the confirmed model directly in R.

------------------------------------------------------------------------

## Full worked example: AVC death after repair

This mirrors the final multivariable model from
`examples/hm.death.AVC.sas` in the reference C repository.

### SAS (original)

``` sas
%HAZARD(
PROC HAZARD DATA=AVCS P CONSERVE OUTHAZ=OUTEST CONDITION=14 QUASI;
     TIME INT_DEAD;
     EVENT DEAD;
     PARMS MUE=0.3504743 THALF=0.1905077 NU=1.437416 M=1 FIXM
           MUC=4.391673E-07;
     EARLY   AGE=-0.03205774, COM_IV=1.336675,  MAL=0.6872028,
             OPMOS=-0.01963377, OP_AGE=0.0002086689, STATUS=0.5169533;
     CONSTANT INC_SURG=1.375285, ORIFICE=3.11765, STATUS=1.054988;
);
```

### R equivalent (current runnable translation pattern)

``` r
# Assumed: avcs is a data.frame read from the AVC flat file
# (same variables as the SAS DATA step)

avcs <- avcs |>
  transform(
    LN_AGE   = log(AGE),
    LN_OPMOS = log(OPMOS),
    LN_INC   = ifelse(is.na(INC_SURG), NA, log(INC_SURG + 1)),
    LN_NYHA  = log(STATUS)
  )

# Replace missing INC_SURG with column mean (mirrors PROC STANDARD REPLACE)
avcs$INC_SURG[is.na(avcs$INC_SURG)] <- mean(avcs$INC_SURG, na.rm = TRUE)

X <- data.matrix(avcs[, c("AGE", "COM_IV", "MAL", "OPMOS", "OP_AGE",
                           "STATUS", "INC_SURG", "ORIFICE")])

fit <- hazard(
  time    = avcs$INT_DEAD,
  status  = avcs$DEAD,
  x       = X,
  theta   = c(
    # Hazard shape parameters (from PARMS)
    MUE   = 0.3504743,
    THALF = 0.1905077,
    NU    = 1.437416,
    M     = 1,
    MUC   = 4.391673e-07,
    # Covariate coefficients (from EARLY + CONSTANT blocks)
    AGE      = -0.03205774,
    COM_IV   =  1.336675,
    MAL      =  0.6872028,
    OPMOS    = -0.01963377,
    OP_AGE   =  0.0002086689,
    STATUS   =  0.5169533,
    INC_SURG =  1.375285,
    ORIFICE  =  3.11765
  ),
  dist    = "weibull",
  control = list(
    condition = 14,
    quasi     = TRUE,
    conserve  = TRUE,
    fix       = c("M")   # FIXM
  )
)

fit
```

------------------------------------------------------------------------

## Data preparation differences

| SAS HAZARD                          | TemporalHazard R                                                                                                                           |
|-------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| `PROC STANDARD REPLACE` for missing | Replace with `mean(..., na.rm = TRUE)` manually                                                                                            |
| Log transforms in `DATA` step       | [`transform()`](https://rdrr.io/r/base/transform.html) or [`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)           |
| Variable labels via `LABEL`         | Column names of `x` matrix carry through                                                                                                   |
| `FILENAME` / `INFILE` / `INPUT`     | [`read.table()`](https://rdrr.io/r/utils/read.table.html), [`read.csv()`](https://rdrr.io/r/utils/read.table.html), or `haven::read_sas()` |
| SAS formats (date, currency, etc.)  | Standard R numeric; apply [`as.Date()`](https://rdrr.io/r/base/as.Date.html) if needed                                                     |

------------------------------------------------------------------------

## Output object

[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
returns a list of class `hazard`:

| Slot             | Contents                                                                                 |
|------------------|------------------------------------------------------------------------------------------|
| `$call`          | Unevaluated [`match.call()`](https://rdrr.io/r/base/match.call.html) for reproducibility |
| `$spec$dist`     | Baseline distribution label                                                              |
| `$spec$control`  | Control options passed at construction                                                   |
| `$data$time`     | Follow-up time vector                                                                    |
| `$data$status`   | Event indicator (numeric)                                                                |
| `$data$x`        | Design matrix                                                                            |
| `$fit$theta`     | Coefficient vector (starting values at M1; fitted at M2+)                                |
| `$fit$converged` | `NA` at M1; `TRUE`/`FALSE` from M2 optimizer                                             |
| `$fit$objective` | Log-likelihood at convergence (M2+)                                                      |
| `$legacy_args`   | Named pass-through arguments for parity                                                  |

------------------------------------------------------------------------

## Prediction

[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
currently supports two output types:

``` r
# Linear predictor (X %*% theta)
eta <- predict(fit, type = "linear_predictor")

# Hazard scale (exp(eta))
hz  <- predict(fit, type = "hazard")
```

From M2 onward, `"survival"` and `"cumulative_hazard"` types will be
added to mirror the `P` (predict/print) option in `PROC HAZARD`.

------------------------------------------------------------------------

## Rcpp acceleration note

`TemporalHazard` is currently pure R. If profiling against large real
datasets reveals bottlenecks in the likelihood kernel or optimizer inner
loop, those specific functions will be re-implemented with Rcpp. The
public interface
([`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md),
[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md),
etc.) will not change.

## References

Blackstone EH, Naftel DC, Turner ME Jr. The decomposition of
time-varying hazard into phases, each incorporating a separate stream of
concomitant information. *J Am Stat Assoc.* 1986;81(395):615–624. doi:
[10.1080/01621459.1986.10478314](https://doi.org/10.1080/01621459.1986.10478314)
