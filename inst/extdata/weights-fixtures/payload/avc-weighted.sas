/* =======================================================================
 * avc-weighted.sas  --  TemporalHazard fractional-weight parity capture
 *
 * Run via run.sh (which sets HAZAPPS + HZ_MACROS before calling SAS):
 *
 *   cd /path/to/this/folder
 *   ./run.sh
 *
 * Produces ./capture/avc-weighted.lst, which parse-weighted-lst.R parses
 * into inst/fixtures/weighted-avc-fractional.rds.
 *
 * Purpose: confirm SAS PROC HAZARD's WEIGHT statement matches R's
 * hazard(weights = ...) for NON-INTEGER (inverse-probability-style)
 * weights -- roadmap 7a / FIXTURE-GAP-LIST B5.  R-side correctness is
 * already proven by the additive-split + linear-scaling invariants in
 * test-weights.R; this fixture is the external SAS confirmation.
 *
 * Model: two-phase (early CDF + constant) hazard with Weibull shaping,
 * AGE as a covariate in both phases, weighted by the fractional column
 * IPW.  Shape parameters (THALF, NU, M) are fixed at the AVC
 * unconditional fit values (same as the stepwise capture) for numerical
 * stability.  Starting MU values likewise.
 * ===================================================================== */

/* -- Load HAZARD macros ------------------------------------------------- */
%let _macros = %sysget(HZ_MACROS);
%if %length(&_macros) = 0 %then %do;
  %put ERROR: HZ_MACROS env var not set -- run this script via run.sh.;
  %abort;
%end;
OPTIONS SASAUTOS=("&_macros" SASAUTOS);

/* Override TMPDIR so %HAZARD macro and C binary agree on the XPORT path. */
%let _hz_tmp = %sysget(HZ_TMPDIR);
%if %length(&_hz_tmp) = 0 %then %do;
  %let _hz_tmp = /tmp;
%end;
OPTIONS SET=TMPDIR "&_hz_tmp";

/* SAS site config overrides HAZAPPS; override it back to the local build. */
%let _hz_apps = %sysget(HAZAPPS);
%if %length(&_hz_apps) > 0 %then %do;
  OPTIONS SET=HAZAPPS "&_hz_apps";
%end;

/* -- Paths ---------------------------------------------------------------- */
%LET csv_path = ./avc-weighted.csv;
%LET out_dir  = ./capture;
%LET run_date = 2026-06-11;

/* -- Import avc-weighted.csv ---------------------------------------------- */
PROC IMPORT DATAFILE = "&csv_path"
            OUT      = _avc
            DBMS     = CSV
            REPLACE;
  GETNAMES = YES;
RUN;

/* Exclude rows with any missing analysis variable (matches R na.omit).
   IPW is a function of AGE, so it is non-missing wherever AGE is.          */
DATA _avc_clean;
  SET _avc;
  IF NMISS(INT_DEAD, DEAD, AGE, IPW) = 0;
RUN;

/* -- Weighted two-phase fit via PROC HAZARD ------------------------------- */
/*                                                                           */
/* WEIGHT IPW applies the fractional column as Fisher observation weights.   */
/* No SELECTION -- AGE is forced into both phases so the model is fixed and  */
/* the R parity refit uses the identical specification.                      */

%let hzrdelds = 1;
%HAZARD(
PROC HAZARD DATA = _avc_clean NOCOV NOCOR CONDITION=14;
     TIME INT_DEAD;
     EVENT DEAD;
     WEIGHT IPW;
     PARMS MUE=0.2361727 THALF=0.1512095 NU=1.438652 M=1 FIXM
           FIXTHALF FIXNU
           MUC=0.0005436977;
     EARLY
       AGE;
     CONSTANT
       AGE;
);

/* -- Metadata file -------------------------------------------------------- */
OPTIONS DLCREATEDIR;
LIBNAME _cap "&out_dir";
LIBNAME _cap CLEAR;

DATA _NULL_;
  FILE "&out_dir/weighted_meta.txt";
  PUT "sas_version=9.4";
  PUT "dataset=avc";
  PUT "dist=multiphase";
  PUT "model=early_cdf+constant";
  PUT "covariates=age";
  PUT "weight_var=ipw";
  PUT "weight_kind=fractional";
  PUT "conserve=1";
  PUT "proc_hazard=SAS HAZARD (PROC HAZARD via %HAZARD macro), WEIGHT IPW";
  PUT "captured_on=&run_date";
RUN;
