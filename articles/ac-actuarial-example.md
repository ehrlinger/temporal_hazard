# Example: Actuarial Estimates (ac.\*)

``` r
library(TemporalHazard)
library(hvtiPlotR)
```

These examples correspond to the `ac.*` SAS files: nonparametric
actuarial (Kaplan-Meier style) event-time estimates used to characterize
the raw data before parametric modeling.

------------------------------------------------------------------------

## ac.death.AVC — Actuarial estimate, AVC death

**SAS source:** `examples/ac.death.AVC.sas`

**Purpose:** Compute the Kaplan-Meier survival estimate for death after
AVC repair. This is the nonparametric foundation before fitting the
parametric hazard model.

### SAS (original)

``` sas
PROC LIFETEST DATA=AVCS METHOD=KM PLOTS=(S H) OUTSURV=KMOUT;
  TIME INT_DEAD * DEAD(0);
RUN;
```

### R translation

``` r
# avcs <- read.table("data/avc", ...)

# Kaplan-Meier using survival package
library(survival)

km_fit <- survfit(Surv(INT_DEAD, DEAD) ~ 1, data = avcs)
summary(km_fit)

# Plot with hvtiPlotR
hv_obj <- hv_survival(avcs, time_col = "INT_DEAD", event_col = "DEAD")
plot(hv_obj, type = "survival") +
  ggplot2::labs(
    x = "Years after repair",
    y = "Freedom from death (%)",
    title = "Kaplan-Meier: Death after AVC repair"
  ) +
  hv_theme()
```

### Notes

- The `ac.*` examples establish the baseline nonparametric estimate
  against which the fitted parametric hazard model is compared visually.
- In HAZARD workflow, you always run `ac.*` first to characterize the
  raw data, then `hz.*` to fit the parametric shape.
- [`hvtiPlotR::hv_survival()`](https://ehrlinger.github.io/hvtiPlotR/reference/hv_survival.html) +
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) produces
  publication-quality KM plots with confidence intervals matching the
  SAS `PROC LIFETEST PLOTS=(S)` output style.
