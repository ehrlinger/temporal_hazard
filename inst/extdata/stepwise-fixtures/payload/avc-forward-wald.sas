/* =======================================================================
 * avc-forward-wald.sas  --  TemporalHazard stepwise parity capture
 *
 * Run via run.sh (which sets HAZAPPS + HZ_MACROS before calling SAS):
 *
 *   cd /path/to/this/folder
 *   ./run.sh
 *
 * Produces ./capture/avc-forward-wald.lst, which the R fixture builder
 * parses to create inst/fixtures/stepwise-avc-forward-wald.rds.
 *
 * Syntax based on hm.death.AVC.sas from the HAZARD example corpus.
 * Starting parameters from the unconditional hz.death.AVC fit.
 * Shaping parameters (THALF, NU, M) are fixed during stepwise for speed
 * and numerical stability (standard practice).
 * ===================================================================== */

/* -- Load HAZARD macros ------------------------------------------------- */
/* HZ_MACROS env var is set by run.sh to the dir containing hazard.sas     */
%let _macros = %sysget(HZ_MACROS);
%if %length(&_macros) = 0 %then %do;
  %put ERROR: HZ_MACROS env var not set -- run this script via run.sh.;
  %abort;
%end;
OPTIONS SASAUTOS=("&_macros" SASAUTOS);

/* Override TMPDIR so %HAZARD macro and C binary agree on the XPORT path.
   SAS Foundation overwrites the shell's $TMPDIR; OPTIONS SET= patches it
   back before %HAZARD writes the handoff file.                            */
%let _hz_tmp = %sysget(HZ_TMPDIR);
%if %length(&_hz_tmp) = 0 %then %do;
  %let _hz_tmp = /tmp;
%end;
OPTIONS SET=TMPDIR "&_hz_tmp";

/* SAS site config overrides HAZAPPS; override it back to the local build.  */
%let _hz_apps = %sysget(HAZAPPS);
%if %length(&_hz_apps) > 0 %then %do;
  OPTIONS SET=HAZAPPS "&_hz_apps";
%end;

/* -- Paths ---------------------------------------------------------------- */
%LET csv_path = ./avc.csv;
%LET out_dir  = ./capture;
%LET run_date = 2026-06-10;

/* -- Import avc.csv ------------------------------------------------------- */
PROC IMPORT DATAFILE = "&csv_path"
            OUT      = _avc
            DBMS     = CSV
            REPLACE;
  GETNAMES = YES;
RUN;

/* Exclude rows with any missing analysis variable (matches R na.omit) */
DATA _avc_clean;
  SET _avc;
  IF NMISS(INT_DEAD, DEAD, AGE, MAL, INC_SURG, ORIFICE, COM_IV) = 0;
RUN;

/* -- Forward Wald stepwise via PROC HAZARD -------------------------------- */
/*                                                                           */
/* Model: two-phase (early + constant) hazard with Weibull shaping.         */
/* Candidates: AGE, MAL, INC_SURG, ORIFICE, COM_IV in both phases.         */
/* SLE=0.30, SLS=0.20 match the R hzr_stepwise() call in the parity test.  */
/*                                                                           */
/* Output goes to the SAS listing (avc-forward-wald.lst); run.sh copies it  */
/* to ./capture/ for the R fixture builder.                                 */

%let hzrdelds = 1;
%HAZARD(
PROC HAZARD DATA = _avc_clean NOCOV NOCOR CONDITION=14;
     TIME INT_DEAD;
     EVENT DEAD;
     PARMS MUE=0.2361727 THALF=0.1512095 NU=1.438652 M=1 FIXM
           FIXTHALF FIXNU
           MUC=0.0005436977;
     SELECTION SLE=0.30 SLS=0.20 MAXSTEPS=50;
     EARLY
       AGE, MAL, INC_SURG, ORIFICE, COM_IV;
     CONSTANT
       AGE, MAL, INC_SURG, ORIFICE, COM_IV;
);

/* -- Metadata file -------------------------------------------------------- */
OPTIONS DLCREATEDIR;
LIBNAME _cap "&out_dir";
LIBNAME _cap CLEAR;

DATA _NULL_;
  FILE "&out_dir/stepwise_meta.txt";
  PUT "sas_version=9.4";
  PUT "dataset=avc";
  PUT "dist=weibull";
  PUT "criterion=wald";
  PUT "direction=forward";
  PUT "slentry=0.30";
  PUT "slstay=0.20";
  PUT "proc_hazard=SAS HAZARD (PROC HAZARD via %HAZARD macro)";
  PUT "captured_on=&run_date";
RUN;
