## ----setup, include=FALSE-----------------------------------------------------
has_pkg <- requireNamespace("TemporalHazard", quietly = TRUE)
if (has_pkg) library(TemporalHazard)
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>",
                      eval = has_pkg)


## ----file-layout, echo=FALSE--------------------------------------------------
layout <- data.frame(
  Layer = c(
    rep("User API", 2),
    rep("Multiphase engine", 3),
    rep("Single-phase likelihoods", 4),
    rep("Shared infrastructure", 5)
  ),
  File = c(
    "hazard_api.R", "argument_mapping.R",
    "likelihood-multiphase.R", "phase-spec.R", "decomposition.R",
    "likelihood-weibull.R", "likelihood-exponential.R",
    "likelihood-loglogistic.R", "likelihood-lognormal.R",
    "optimizer.R", "math_primitives.R", "formula-helpers.R",
    "golden_fixtures.R", "parity-helpers.R"
  ),
  Purpose = c(
    "hazard(), predict.hazard(), print/summary/coef/vcov S3 methods",
    "SAS HAZARD -> R parameter translation table",
    "Additive N-phase likelihood, cumhaz/hazard evaluators, multi-start optimizer",
    "hzr_phase() constructor, validators, starting-value assembly",
    "Unified decompos(t; t_half, nu, m) parametric family",
    "Weibull PH likelihood, gradient, optimizer",
    "Exponential (constant hazard) likelihood",
    "Log-logistic (proportional odds) likelihood",
    "Log-normal (AFT) likelihood",
    "Generic L-BFGS-B/BFGS optimizer with Hessian-based vcov",
    "Numerically stable log1pexp, log1mexp, clamp_prob",
    "Surv() formula parsing for right/left/interval censoring",
    "Synthetic fixture generators (.rds reference outputs)",
    "Stubs for cross-validating against C HAZARD binary"
  )
)
knitr::kable(layout, caption = "R source files by architectural layer")


## ----theta-layout, echo=FALSE-------------------------------------------------
theta_layout <- data.frame(
  `Phase type` = c("cdf / hazard", "g3", "constant"),
  `Sub-vector` = c(
    "[log_mu, log_t_half, nu, m, beta_1, ..., beta_p]",
    "[log_mu, log_tau, gamma, alpha, eta, beta_1, ..., beta_p]",
    "[log_mu, beta_1, ..., beta_p]"
  ),
  `Count` = c("4 + p", "5 + p", "1 + p")
)
knitr::kable(theta_layout, caption = "Internal theta layout per phase")


## ----fixtures, echo=FALSE-----------------------------------------------------
fixtures <- data.frame(
  File = c(
    "hz_univariate.rds",
    "hm_multivariate.rds",
    "hm_edge_case.rds",
    "hz_loglogistic.rds",
    "hz_lognormal.rds",
    "mp_synthetic_3phase.rds",
    "mp_c_reference_kul.rds"
  ),
  Distribution = c(
    "Weibull", "Weibull", "Weibull",
    "Log-logistic", "Log-normal",
    "Multiphase (3)", "Multiphase (3)"
  ),
  Description = c(
    "Univariable shape estimation (n=100)",
    "2 covariates (n=100)",
    "Edge case: n=20, 3 covariates",
    "Univariable (n=80)",
    "Univariable (n=80)",
    "Synthetic early CDF + constant + late hazard (n=500)",
    "KUL CABG dataset + C HAZARD binary reference output (n=5880)"
  ),
  Source = c(
    rep("Synthetic (seed=42)", 5),
    "Synthetic (seed=42)",
    "Clinical + C binary"
  )
)
knitr::kable(fixtures, caption = "Golden fixture inventory")


## ----test-tiers, echo=FALSE---------------------------------------------------
tiers <- data.frame(
  Tier = c(
    "Unit tests",
    "Distribution tests",
    "Integration tests",
    "Parity tests"
  ),
  Files = c(
    "test-math-primitives, test-decomposition, test-phase-spec, test-argument-mapping",
    paste("test-gradient-weibull, test-exponential-dist,",
          "test-loglogistic-dist, test-lognormal-dist, test-multiphase-gradient"),
    paste("test-hazard-api, test-predict-types, test-interval-censoring-*,",
          "test-time-varying-*, test-multiphase-likelihood, test-multiphase-api"),
    "test-parity-core, test-parity-edge-cases, test-parity-c-binary, test-multiphase-parity"
  ),
  Purpose = c(
    "Verify individual functions in isolation",
    "Likelihood, gradient, and optimizer for each distribution",
    "End-to-end: hazard() -> predict() -> summary() pipeline, censoring types, multiphase wiring",
    "Golden fixture round-trip, C binary cross-validation"
  )
)
knitr::kable(tiers, caption = "Test suite tiers")


## ----datasets, echo=FALSE-----------------------------------------------------
datasets <- data.frame(
  File = c("avc.csv", "cabgkul.csv", "omc.csv", "tga.csv", "valves.csv"),
  Study = c(
    "Atrioventricular canal repair",
    "Coronary artery bypass grafting (KU Leuven)",
    "Open mitral commissurotomy",
    "Transposition of great arteries (arterial switch)",
    "Primary valve replacement"
  ),
  n = c(310, 5880, 339, 470, 1533),
  Events = c("70 deaths", "545 deaths", "thromboembolic events (repeated)", "deaths", "deaths, PVE, reoperation"),
  Covariates = c(
    "NYHA, age, anatomy, era",
    "None (intercept only)",
    "TE events (repeated measures)",
    "anatomy, coronary pattern, era",
    "age, NYHA, valve position, pathology, race"
  ),
  `SAS origin` = c(
    "hz.death.AVC, hm.death.AVC",
    "hz.deadp.KUL",
    "hz.te123.OMC",
    "hs.dthar.TGA",
    "hm.deadp.VALVES"
  )
)
knitr::kable(datasets, caption = "Reference datasets (lazy-loaded via `data()`)")


## ----load-demo----------------------------------------------------------------
# Datasets are lazy-loaded with the package — just reference them directly.
# Raw CSVs are also available in inst/extdata/ for advanced use:
#   read.csv(system.file("extdata", "cabgkul.csv", package = "TemporalHazard"))

data(cabgkul)
str(cabgkul)

data(avc)
str(avc)


## ----avc-detail, echo=FALSE---------------------------------------------------
avc_vars <- data.frame(
  Variable = c("study", "status", "inc_surg", "opmos", "age", "mal",
               "com_iv", "orifice", "dead", "int_dead", "op_age"),
  Label = c("Study number", "NYHA functional class (I-V)", "Surgical grade of AV valve incompetence",
            "Date of operation (months since Jan 1967)", "Age (months) at repair",
            "Important associated cardiac anomaly", "Interventricular communication",
            "Accessory left AV valve orifice", "Death (event indicator)",
            "Follow-up interval (months)", "Interaction: opmos x age"),
  Type = c("character", rep("numeric", 10))
)
knitr::kable(avc_vars, caption = "AVC dataset variables")


## ----mapping------------------------------------------------------------------
mapping <- hzr_argument_mapping()
knitr::kable(
  mapping[, c("legacy_input", "r_parameter", "implementation_status", "notes")],
  caption = "SAS HAZARD to R parameter mapping (excerpt)"
)

