#' AVC: Atrioventricular Canal Repair
#'
#' Survival data for 310 patients who underwent repair of atrioventricular
#' septal defects (congenital heart disease) at the Cleveland Clinic between
#' 1977 and 1993. Exhibits two identifiable hazard phases: an early
#' post-operative risk and a constant late phase.
#'
#' @format A data frame with 310 rows and 11 variables:
#' \describe{
#'   \item{study}{Patient identifier}
#'   \item{status}{NYHA functional class (1--4)}
#'   \item{inc_surg}{Surgical grade of AV valve incompetence}
#'   \item{opmos}{Date of operation (months since January 1967)}
#'   \item{age}{Age at repair (months)}
#'   \item{mal}{Malalignment indicator (0/1)}
#'   \item{com_iv}{Interventricular communication indicator (0/1)}
#'   \item{orifice}{Associated cardiac anomaly indicator (0/1)}
#'   \item{dead}{Death indicator (1 = dead, 0 = censored)}
#'   \item{int_dead}{Follow-up interval to death or last contact (months)}
#'   \item{op_age}{Interaction term: opmos x age}
#' }
#'
#' @source Blackstone, Naftel, and Turner (1986)
#'   \doi{10.1080/01621459.1986.10478314}. Cleveland Clinic Foundation.
#' @family datasets
"avc"

#' CABGKUL: Primary Isolated Coronary Artery Bypass Grafting (KU Leuven)
#'
#' Survival data for 5,880 patients who underwent primary isolated CABG
#' at KU Leuven, Belgium, between 1971 and July 1987. The simplest dataset
#' structure (intercept-only, right-censored) with large sample size
#' exercising all three temporal hazard phases.
#'
#' @format A data frame with 5880 rows and 2 variables:
#' \describe{
#'   \item{int_dead}{Follow-up interval to death or last contact (months)}
#'   \item{dead}{Death indicator (1 = dead, 0 = censored)}
#' }
#'
#' @source KU Leuven cardiac surgery registry. Primary benchmark dataset for
#'   C binary parity testing.
#' @family datasets
"cabgkul"

#' OMC: Open Mitral Commissurotomy
#'
#' Data for 339 patients who underwent open mitral commissurotomy at the
#' University of Alabama Birmingham. Contains repeated thromboembolic events
#' (up to 3 per patient) with left censoring, exercising the interval
#' censoring likelihood.
#'
#' @format A data frame with 339 rows and 7 variables:
#' \describe{
#'   \item{study}{Patient identifier}
#'   \item{te1}{Indicator for first thromboembolic event}
#'   \item{te2}{Indicator for second thromboembolic event}
#'   \item{te3}{Indicator for third thromboembolic event}
#'   \item{int_dead}{Follow-up interval to death or last contact (months)}
#'   \item{dead}{Death indicator (1 = dead, 0 = censored)}
#'   \item{opdjul}{Operation date (Julian)}
#' }
#'
#' @source University of Alabama Birmingham cardiac surgery registry.
#' @family datasets
"omc"

#' TGA: Transposition of the Great Arteries
#'
#' Survival data for 470 patients who underwent the arterial switch operation
#' for transposition of the great arteries at Boston Children's Hospital and
#' Children's Hospital of Philadelphia. Used for sensitivity analysis and
#' internal validation examples.
#'
#' @format A data frame with 470 rows and 14 variables:
#' \describe{
#'   \item{study}{Patient identifier}
#'   \item{simple}{Simple TGA indicator (0/1)}
#'   \item{dextroin}{D-looped transposition indicator (0/1)}
#'   \item{ca_1rl2c}{Coronary artery pattern indicator}
#'   \item{hyaaproc}{Hybrid approach procedure indicator (0/1)}
#'   \item{no_tca}{No total circulatory arrest indicator (0/1)}
#'   \item{tca_time}{Total circulatory arrest time (minutes)}
#'   \item{age_days}{Age at operation (days)}
#'   \item{arciopyr}{Aortic cross-clamp time per year}
#'   \item{dead}{Death indicator (1 = dead, 0 = censored)}
#'   \item{int_dead}{Follow-up interval to death or last contact (months)}
#'   \item{source}{Source institution (BCH or CHOP)}
#'   \item{ca1_2_l}{Coronary artery configuration (1/2/L)}
#'   \item{opyear}{Year of operation}
#' }
#'
#' @source Boston Children's Hospital and Children's Hospital of Philadelphia.
#' @family datasets
"tga"

#' Valves: Primary Heart Valve Replacement
#'
#' Data for 1,533 patients who underwent primary heart valve replacement.
#' The largest multivariable example dataset with multiple endpoints
#' including death, prosthetic valve endocarditis (PVE), bioprosthesis
#' degeneration, and reoperation.
#'
#' @format A data frame with 1533 rows and 19 variables:
#' \describe{
#'   \item{age_cop}{Age at operation (years)}
#'   \item{nyha}{NYHA functional class (1--4)}
#'   \item{mitral}{Mitral valve position indicator (0/1)}
#'   \item{double_}{Double valve replacement indicator (0/1)}
#'   \item{ao_pinc}{Aortic position, incompetence (0/1)}
#'   \item{black}{Black race indicator (0/1)}
#'   \item{i_path}{Ischemic pathology indicator (0/1)}
#'   \item{nve}{Native valve endocarditis indicator (0/1)}
#'   \item{mechvalv}{Mechanical valve indicator (0/1)}
#'   \item{male}{Male sex indicator (0/1)}
#'   \item{int_dead}{Follow-up interval to death or last contact (months)}
#'   \item{dead}{Death indicator (1 = dead, 0 = censored)}
#'   \item{int_pve}{Follow-up interval to PVE or last contact (months)}
#'   \item{pve}{PVE indicator (1 = PVE, 0 = censored)}
#'   \item{bio}{Bioprosthesis indicator (0/1)}
#'   \item{int_rdg}{Follow-up interval to degeneration or last contact (months)}
#'   \item{reop_dg}{Reoperation for degeneration indicator (0/1)}
#'   \item{int_reop}{Follow-up interval to reoperation or last contact (months)}
#'   \item{reop}{Reoperation indicator (0/1)}
#' }
#'
#' @source Cleveland Clinic Foundation heart valve replacement registry.
#' @family datasets
"valves"
