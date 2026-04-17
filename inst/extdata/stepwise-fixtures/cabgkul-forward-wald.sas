/* =======================================================================
 * cabgkul-forward-wald.sas — Reference capture for R-side stepwise parity
 *
 * Runs forward Wald stepwise selection on CABGKUL using SAS PROC HAZARD
 * with SLENTRY = 0.30 / SLSTAY = 0.20 (SAS defaults) and emits three
 * text files consumed by .hzr_build_stepwise_fixture() in R.
 *
 * Usage:
 *   1. Point LIBNAME to the directory containing cabgkul.sas7bdat.
 *   2. Adjust %LET out_dir if you want the capture written elsewhere.
 *   3. Run in SAS 9.4+.  The three output files can then be copied
 *      into the R session for fixture construction.
 * ===================================================================== */

%LET in_dir  = /path/to/cabgkul;
%LET out_dir = /path/to/capture;

LIBNAME hvti "&in_dir";

/* -- PROC HAZARD invocation --------------------------------------------- */
ODS LISTING CLOSE;
ODS OUTPUT
  SelectionSummary = _sel_trace     /* step-by-step history         */
  ParameterEstimates = _sel_final   /* final coefficient table      */
  ConvergenceStatus = _sel_conv     /* iterations, convergence code */
  FitStatistics = _sel_fit;         /* logLik, AIC, BIC             */

PROC HAZARD DATA = hvti.cabgkul
            METHOD = STEPWISE
            SLENTRY = 0.30
            SLSTAY  = 0.20;
  TIME INT_DEAD * DEAD(0);
  MODEL INT_DEAD = AGE MALE INC_SURG ORIFICE COM_IV
        / DIST = WEIBULL;
RUN;

ODS LISTING;

/* -- Write the three capture files -------------------------------------- */

/* 1. Step trace: one row per step. */
DATA _step_trace;
  SET _sel_trace;
  LENGTH action $ 8 variable $ 32 phase $ 8;
  step_num = _N_;
  /* Map SAS "Effect Entered" / "Effect Removed" columns → action */
  IF MISSING(EffectRemoved) THEN DO;
    action   = "enter";
    variable = EffectEntered;
  END;
  ELSE DO;
    action   = "drop";
    variable = EffectRemoved;
  END;
  phase    = ".";                      /* placeholder for single-dist */
  stat     = ChiSquare;
  df       = DF;
  p_value  = ProbChiSq;
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
  variable  = Parameter;
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

/* 3. Meta / GOF totals — a tiny key=value text file. */
DATA _NULL_;
  FILE "&out_dir/stepwise_meta.txt";
  SET _sel_fit(WHERE = (UPCASE(Criterion) = "-2 LOG L" OR
                         UPCASE(Criterion) = "AIC"));
  IF UPCASE(Criterion) = "-2 LOG L" THEN DO;
    logLik_twice = Value;
    PUT "logLik=" (logLik_twice / -2) 12.6;
  END;
  ELSE IF UPCASE(Criterion) = "AIC" THEN DO;
    PUT "aic=" Value 12.6;
  END;
RUN;

DATA _NULL_;
  FILE "&out_dir/stepwise_meta.txt" MOD;
  SET _sel_conv;
  PUT "iterations=" _N_;
  PUT "converged=" Status;
  PUT "sas_version=9.4";
  PUT "dataset=cabgkul";
  PUT "dist=weibull";
  PUT "criterion=wald";
  PUT "direction=forward";
  PUT "slentry=0.30";
  PUT "slstay=0.20";
  PUT "captured_on=2026-04-17";
RUN;
