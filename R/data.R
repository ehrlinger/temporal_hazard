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
#'
#' @examples
#' data(avc)
#' avc <- na.omit(avc)
#'
#' # Kaplan-Meier survival
#' km <- survival::survfit(survival::Surv(int_dead, dead) ~ 1, data = avc)
#' plot(km, xlab = "Months after AVC repair", ylab = "Survival",
#'      main = "AVC: Kaplan-Meier survival estimate")
#'
#' \donttest{
#' # Multiphase hazard fit
#' fit <- hazard(
#'   survival::Surv(int_dead, dead) ~ 1, data = avc,
#'   dist = "multiphase",
#'   phases = list(
#'     early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1),
#'     constant = hzr_phase("constant"),
#'     late     = hzr_phase("cdf", t_half = 10, nu = 1, m = 1)
#'   ),
#'   fit = TRUE, control = list(n_starts = 5, maxit = 1000)
#' )
#' summary(fit)
#' }
#'
#' @seealso \code{vignette("fitting-hazard-models")},
#'   \code{vignette("prediction-visualization")}
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
#'
#' @examples
#' data(cabgkul)
#'
#' # Kaplan-Meier survival
#' km <- survival::survfit(survival::Surv(int_dead, dead) ~ 1, data = cabgkul)
#' plot(km, xlab = "Months after CABG", ylab = "Survival",
#'      main = "CABGKUL: Kaplan-Meier survival (n = 5,880)")
#'
#' \donttest{
#' # Single-phase Weibull fit with parametric overlay
#' fit <- hazard(survival::Surv(int_dead, dead) ~ 1, data = cabgkul,
#'               dist = "weibull", theta = c(mu = 0.10, nu = 1.0), fit = TRUE)
#' t_grid <- seq(0.01, max(cabgkul$int_dead) * 0.9, length.out = 200)
#' surv   <- predict(fit, newdata = data.frame(time = t_grid),
#'                   type = "survival")
#' plot(km, xlab = "Months after CABG", ylab = "Survival",
#'      main = "CABGKUL: Weibull vs. Kaplan-Meier")
#' lines(t_grid, surv, col = "blue", lwd = 2)
#' legend("bottomleft", c("KM", "Weibull"), col = c("black", "blue"),
#'        lty = 1, lwd = c(1, 2))
#' }
#'
#' @seealso \code{vignette("fitting-hazard-models")}
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
#'
#' @examples
#' data(valves)
#' valves_cc <- na.omit(valves)
#'
#' # Kaplan-Meier for two endpoints
#' km_death <- survival::survfit(
#'   survival::Surv(int_dead, dead) ~ 1, data = valves_cc)
#' km_pve <- survival::survfit(
#'   survival::Surv(int_pve, pve) ~ 1, data = valves_cc)
#'
#' plot(km_death, xlab = "Months after valve replacement", ylab = "Survival",
#'      main = "Valves: Death and PVE endpoints")
#' lines(km_pve, col = "red")
#' legend("bottomleft", c("Death", "PVE"), col = c("black", "red"), lty = 1)
#'
#' @seealso \code{vignette("fitting-hazard-models")},
#'   \code{vignette("prediction-visualization")}
#' @family datasets
"valves"
