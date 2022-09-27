#' Simulated example dataset of SNPs obtained during a GWAS
#'
#' A dataset containing the numeric interpretation of 500 SNPs from 200 simulated humans. Both
#' dimensions have a number-based names attribute.
#'
#' @format A matrix with 200 rows and 300 variables with integer values being either 0, 1 or 2.
"example_snp_data"


#' Simulated example dataset for a survival GWAS
#'
#' A dataset containing time-to-event outcomes for 200 subjects, simulated from a proportional
#' hazards model with one interaction increasing the hazard. The censoring times have a uniform
#' distribution.
#'
#' @format A Surv-class object with 200 rows and 2 variables:
#' \describe{
#'   \item{time}{the time until the event took place}
#'   \item{status}{indicator whether or not the subject was censored. 1 if the event was observed.}
#' }
"example_survival_data"
