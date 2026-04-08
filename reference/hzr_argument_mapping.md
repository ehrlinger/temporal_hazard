# Legacy HAZARD to hvtiRhazard argument mapping

Returns a formal mapping table that defines how legacy SAS
HAZARD/C-style inputs map to `hazard(...)` arguments in this package.

## Usage

``` r
hzr_argument_mapping(include_planned = TRUE)
```

## Arguments

- include_planned:

  Logical; if `FALSE`, only rows marked as implemented are returned.

## Value

A data frame with one row per mapping rule.

## Examples

``` r
hzr_argument_mapping()
#>    sas_statement              legacy_input                   c_concept
#> 1         HAZARD             TIME variable              obs time array
#> 2         HAZARD     EVENT/censor variable       event indicator array
#> 3         HAZARD         X covariate block               design matrix
#> 4         HAZARD        initial parameters            parameter vector
#> 5         HAZARD     baseline distribution phase distribution selector
#> 6         HAZARD           control options    optimizer/control struct
#> 7         HAZARD additional legacy options        misc legacy switches
#> 8           TIME                         t                 time vector
#> 9          EVENT                    status                event vector
#> 10         PARMS                    theta0               starting coef
#> 11          DIST                      dist               dist selector
#>    r_parameter required                expected_type
#> 1         time     TRUE               numeric vector
#> 2       status     TRUE       numeric/logical vector
#> 3            x    FALSE numeric matrix or data.frame
#> 4        theta    FALSE               numeric vector
#> 5         dist    FALSE             character scalar
#> 6      control    FALSE                   named list
#> 7          ...    FALSE              named arguments
#> 8         time     TRUE               numeric vector
#> 9       status     TRUE       numeric/logical vector
#> 10       theta    FALSE               numeric vector
#> 11        dist    FALSE             character scalar
#>                                 transform_rule implementation_status
#> 1                      pass through as numeric           implemented
#> 2                        coerce to numeric 0/1           implemented
#> 3                    data.frame -> data.matrix           implemented
#> 4  length must equal ncol(x) when x is present           implemented
#> 5                  normalized lower-case label           implemented
#> 6                       stored in spec$control           implemented
#> 7             stored in legacy_args for parity           implemented
#> 8                                 pass through           implemented
#> 9                            coerce to numeric           implemented
#> 10                  map PARMS/INITIAL to theta               planned
#> 11                           map DIST= to dist           implemented
#>                                                                   notes
#> 1                                          Core observation time input.
#> 2  Event indicator currently retained as numeric in object$data$status.
#> 3          Future versions will support richer design encoding helpers.
#> 4                         Used by predict.hazard as coefficient vector.
#> 5                   Current default is 'weibull'; more options planned.
#> 6             Control list is stored and reserved for optimizer parity.
#> 7          Supports legacy-style pass-through options during migration.
#> 8                           Canonical SAS migration uses TIME= mapping.
#> 9                          Canonical SAS migration uses EVENT= mapping.
#> 10                         SAS PARMS syntax parser not yet implemented.
#> 11                              SAS DIST keyword maps directly to dist.
hzr_argument_mapping(include_planned = FALSE)
#>    sas_statement              legacy_input                   c_concept
#> 1         HAZARD             TIME variable              obs time array
#> 2         HAZARD     EVENT/censor variable       event indicator array
#> 3         HAZARD         X covariate block               design matrix
#> 4         HAZARD        initial parameters            parameter vector
#> 5         HAZARD     baseline distribution phase distribution selector
#> 6         HAZARD           control options    optimizer/control struct
#> 7         HAZARD additional legacy options        misc legacy switches
#> 8           TIME                         t                 time vector
#> 9          EVENT                    status                event vector
#> 10          DIST                      dist               dist selector
#>    r_parameter required                expected_type
#> 1         time     TRUE               numeric vector
#> 2       status     TRUE       numeric/logical vector
#> 3            x    FALSE numeric matrix or data.frame
#> 4        theta    FALSE               numeric vector
#> 5         dist    FALSE             character scalar
#> 6      control    FALSE                   named list
#> 7          ...    FALSE              named arguments
#> 8         time     TRUE               numeric vector
#> 9       status     TRUE       numeric/logical vector
#> 10        dist    FALSE             character scalar
#>                                 transform_rule implementation_status
#> 1                      pass through as numeric           implemented
#> 2                        coerce to numeric 0/1           implemented
#> 3                    data.frame -> data.matrix           implemented
#> 4  length must equal ncol(x) when x is present           implemented
#> 5                  normalized lower-case label           implemented
#> 6                       stored in spec$control           implemented
#> 7             stored in legacy_args for parity           implemented
#> 8                                 pass through           implemented
#> 9                            coerce to numeric           implemented
#> 10                           map DIST= to dist           implemented
#>                                                                   notes
#> 1                                          Core observation time input.
#> 2  Event indicator currently retained as numeric in object$data$status.
#> 3          Future versions will support richer design encoding helpers.
#> 4                         Used by predict.hazard as coefficient vector.
#> 5                   Current default is 'weibull'; more options planned.
#> 6             Control list is stored and reserved for optimizer parity.
#> 7          Supports legacy-style pass-through options during migration.
#> 8                           Canonical SAS migration uses TIME= mapping.
#> 9                          Canonical SAS migration uses EVENT= mapping.
#> 10                              SAS DIST keyword maps directly to dist.
```
