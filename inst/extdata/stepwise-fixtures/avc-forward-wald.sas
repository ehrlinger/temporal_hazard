/* =======================================================================
 * avc-forward-wald.sas — Reference capture for R-side stepwise parity
 *
 * Runs forward Wald stepwise selection on the AVC dataset using SAS
 * PROC HAZARD with SLENTRY = 0.30 / SLSTAY = 0.20 (SAS defaults) and
 * emits three text files consumed by .hzr_build_stepwise_fixture() in R.
 *
 * Dataset: hz.death.AVC  (n = 310, 68 deaths after na.omit)
 * Covariates: AGE MAL INC_SURG ORIFICE COM_IV
 * Distribution: Weibull
 * Direction: Forward
 * Criterion: Wald chi-square
 *
 * Usage:
 *   1. Upload inst/extdata/avc.csv to the server (or point to an
 *      existing SAS dataset).
 *   2. Set the three macro variables below.
 *   3. Run in SAS 9.4+ with PROC HAZARD available.
 *   4. Copy the three output files back to your R session and run:
 *
 *        TemporalHazard:::.hzr_build_stepwise_fixture(
 *          trace_csv = "path/to/stepwise_trace.csv",
 *          final_csv = "path/to/stepwise_final.csv",
 *          meta_txt  = "path/to/stepwise_meta.txt",
 *          out_path  = "inst/fixtures/stepwise-avc-forward-wald.rds"
 *        )
 *
 * ===================================================================== */

%LET csv_path = /path/to/avc.csv;       /* path to uploaded avc.csv  */
%LET out_dir  = /path/to/capture;       /* directory for output files */
%LET run_date = 2026-06-10;             /* update to actual run date  */

/* -- Import the CSV ------------------------------------------------------ */
PROC IMPORT DATAFILE = "&csv_path"
            OUT      = _avc
            DBMS     = CSV
            REPLACE;
  GETNAMES = YES;
RUN;

/* Remove incomplete cases (matches R's na.omit(avc)) */
DATA _avc_clean;
  SET _avc;
  IF NMISS(INT_DEAD, DEAD, AGE, MAL, INC_SURG, ORIFICE, COM_IV) = 0;
RUN;

/* -- PROC HAZARD invocation --------------------------------------------- */
ODS LISTING CLOSE;
ODS OUTPUT
  SelectionSummary  = _sel_trace    /* step-by-step history         */
  ParameterEstimates = _sel_final   /* final coefficient table      */
  ConvergenceStatus  = _sel_conv    /* iterations, convergence code */
  FitStatistics      = _sel_fit;    /* logLik, AIC, BIC             */

PROC HAZARD DATA = _avc_clean
            METHOD = STEPWISE
            SLENTRY = 0.30
            SLSTAY  = 0.20;
  TIME INT_DEAD * DEAD(0);
  MODEL INT_DEAD = AGE MAL INC_SURG ORIFICE COM_IV
        / DIST = WEIBULL;
RUN;

ODS LISTING;

/* -- Write the three capture files -------------------------------------- */

/* 1. Step trace: one row per step. */
DATA _step_trace;
  SET _sel_trace;
  LENGTH action $ 8 variable $ 32 phase $ 8;
  step_num = _N_;
  IF MISSING(EffectRemoved) THEN DO;
    action   = "enter";
    variable = LOWCASE(STRIP(EffectEntered));
  END;
  ELSE DO;
    action   = "drop";
    variable = LOWCASE(STRIP(EffectRemoved));
  END;
  phase   = ".";
  stat    = ChiSquare;
  df      = DF;
  p_value = ProbChiSq;
  KEEP step_num action variable phase stat df p_value;
RUN;

PROC EXPORT DATA = _step_trace
            OUTFILE = "&out_dir/stepwise_trace.csv"
            DBMS = CSV REPLACE;
RUN;

/* 2. Final coefficient table. */
DATA _final_coef;
  SET _sel_final;
  LENGTH variable $ 32 phase $ 8;
  variable  = LOWCASE(STRIP(Parameter));
  phase     = ".";
  estimate  = Estimate;
  std_error = StdErr;
  z_stat    = Estimate / StdErr;
  p_value   = ProbChiSq;
  KEEP variable phase estimate std_error z_stat p_value;
RUN;

PROC EXPORT DATA = _final_coef
            OUTFILE = "&out_dir/stepwise_final.csv"
            DBMS = CSV REPLACE;
RUN;

/* 3. Meta / GOF totals. */
DATA _NULL_;
  FILE "&out_dir/stepwise_meta.txt";
  SET _sel_fit(WHERE = (UPCASE(STRIP(Criterion)) = "-2 LOG L" OR
                        UPCASE(STRIP(Criterion)) = "AIC"));
  IF UPCASE(STRIP(Criterion)) = "-2 LOG L" THEN
    PUT "logLik=" (Value / -2) 12.6;
  ELSE IF UPCASE(STRIP(Criterion)) = "AIC" THEN
    PUT "aic=" Value 12.6;
RUN;

DATA _NULL_;
  FILE "&out_dir/stepwise_meta.txt" MOD;
  SET _sel_conv(OBS = 1);
  PUT "iterations=" _N_;
  PUT "converged=" Status;
RUN;

DATA _NULL_;
  FILE "&out_dir/stepwise_meta.txt" MOD;
  PUT "sas_version=9.4";
  PUT "dataset=avc";
  PUT "dist=weibull";
  PUT "criterion=wald";
  PUT "direction=forward";
  PUT "slentry=0.30";
  PUT "slstay=0.20";
  PUT "proc_hazard=SAS HAZARD";
  PUT "captured_on=&run_date";
RUN;
