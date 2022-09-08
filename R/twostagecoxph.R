#' Perform a two stage analysis on a survival dataset to detect epistasis
#'
#' @param survival.dataset The survival dataset describing the outcome.
#' @param covariate.matrix The nxp-matrix of covariates of the p covariates of the n patients.
#' @param first.stage.threshold numeric scalar denoting the threshold for the first stage. If a covariate
#'          has a p-value lower than this threshold, it will be passed on to the second stage.
#' @param multiple.hypotheses.correction Correction method, a character string. Passed to stats::p.adjust
#' @param multicore Currently unused
#'
#' @return Hopefully a nicely structured list of some kind?
#' @export
#'
#' @note Parallel processing is currently not supported
#'
#' @examples
#' survival.dataset <- survival::Surv(c(5,5,3,3,2,2,2,1,1,1),
#'                                    c(0,0,1,1,1,1,1,1,1,1))
#' covariate.matrix <- matrix(c(2,2,1,
#'                              2,2,1,
#'                              1,2,1,
#'                              2,1,1,
#'                              1,1,1,
#'                              1,1,1,
#'                              1,0,0,
#'                              0,1,0,
#'                              1,0,0,
#'                              0,0,0),
#'                            nrow = 10, ncol = 3, byrow = TRUE)
#' twostagecoxph(survival.dataset, covariate.matrix)
twostagecoxph <- function(survival.dataset, covariate.matrix, first.stage.threshold = 0.05,
                          multiple.hypotheses.correction = "bonferroni", multicore = FALSE){
  if(!survival::is.Surv(survival.dataset)) stop("survival.dataset must be a Surv object")
  if(!is.matrix(covariate.matrix)) stop("covariate.matrix must be a (strict) matrix")

  if(!(dim(survival.dataset)[1] %in% dim(covariate.matrix))) stop("No dimension of covariate.matrix matches dimension of survival.dataset")
  if(dim(survival.dataset)[1] != dim(covariate.matrix)[1] &&
     dim(survival.dataset)[1] == dim(covariate.matrix)[2]) {
    warning("Transposed the covariate matrix to match dimension with survival.dataset")
    covariate.matrix <- t(covariate.matrix)
  }


  if(!is.numeric(first.stage.threshold)) stop("first.stage.threshold must be numeric")
  if(length(first.stage.threshold) != 1) stop("first.stage.threshold must be of length 1")
  if(first.stage.threshold > 1 || first.stage.threshold < 0) stop("first.stage.threshold must be in the interval [0,1]")

  if(multicore != FALSE) stop("Package currently does not support parallel processing")

  if(multicore == FALSE){
    return(singlecore.twostagecoxph(survival.dataset, covariate.matrix,
                                    first.stage.threshold, multiple.hypotheses.correction))
  }


}


singlecore.twostagecoxph <- function(survival.dataset, covariate.matrix, first.stage.threshold, multiple.hypotheses.correction){
  first.stage.result <- firststagecoxph(survival.dataset, covariate.matrix)

  first.stage.rejections <- (first.stage.result < first.stage.threshold)
  passed.indices <- which(first.stage.rejections == TRUE)
  amount.rejections <- sum(first.stage.rejections)

  second.stage.result.matrix <- matrix(1, nrow = amount.rejections, ncol = amount.rejections,
                                       dimnames = list(passed.indices, passed.indices))
  for (first.index in 1:(amount.rejections-1)){
    for (second.index in (first.index+1):amount.rejections){
      index.first.covariate  <- passed.indices[first.index]
      index.second.covariate <- passed.indices[second.index]
      covariate.one <- covariate.matrix[, index.first.covariate]
      covariate.two <- covariate.matrix[, index.second.covariate]

      fitted.model <- survival::coxph(survival.dataset ~ covariate.one * covariate.two)
      second.stage.result.matrix[first.index, second.index] <- summary(fitted.model)$coefficients[3,5]
    }
  }

  corrected.p.values <- stats::p.adjust(second.stage.result.matrix[upper.tri(second.stage.result.matrix)],
                                        method = multiple.hypotheses.correction)
  corrected.result.matrix <- second.stage.result.matrix
  corrected.result.matrix[upper.tri(corrected.result.matrix)] <- corrected.p.values

  return(list(raw.p.value.matrix = second.stage.result.matrix,
              corrected.p.value.matrix = corrected.result.matrix))
}


firststagecoxph <- function(survival.dataset, covariate.matrix){
  p.value.vector <- rep(1, length = dim(covariate.matrix)[2])
  for(covariate.index in 1:length(p.value.vector)){
    this.covariate <- covariate.matrix[, covariate.index]
    fitted.model <- survival::coxph(survival.dataset ~ this.covariate)
    p.value.vector[covariate.index] <- summary(fitted.model)$coefficients[1,5]
  }
  return(p.value.vector)
}


