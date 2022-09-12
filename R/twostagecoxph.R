#' (WIP) Perform a two stage analysis on a survival dataset to detect epistasis
#'
#' Work in progress. Only the most basic functionality works.
#'
#' @param survival.dataset The survival dataset describing the outcome.
#' @param covariate.matrix The nxp-matrix of covariates of the p covariates of the n patients.
#' @param first.stage.threshold numeric scalar denoting the threshold for the first stage. If a covariate
#'          has a p-value lower than this threshold, it will be passed on to the second stage.
#' @param multiple.hypotheses.correction Correction method, a character string. Passed to stats::p.adjust
#' @param multicore Currently unused
#' @param report.lowest.amount integer scalar denoting how many of the most significant interactions
#'                               the function should report
#' @param return.raw logical, default FALSE; whether or not the output should contain the raw p-values or
#'                     the multiple hypotheses corrected p-values.
#' @param progress numeric, default 1000;  how many iterations should pass silently until an update is given
#'                   about runtime and progress until completion of stages. Set to 0 for no output.
#'
#' @return Hopefully a nicely structured list of some kind?
#' @export
#'
#' @note Parallel processing is currently not supported
#'
#' @seealso print.twostageGWAS()
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
#' twostagecoxph(survival.dataset, covariate.matrix, progress = 0)
twostagecoxph <- function(survival.dataset, covariate.matrix, first.stage.threshold = 0.05,
                          multiple.hypotheses.correction = "bonferroni", multicore = FALSE,
                          report.lowest.amount = 5, return.raw = FALSE,
                          progress = 1000){
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
  if(first.stage.threshold == 1) stop("first.stage.threshold is 1, consider decreasing this. The two stage method is not suitable for this.")
  if(first.stage.threshold == 0){
    if(stats::runif(1) < 0.5) {
      stop("first.stage.threshold is 0, consider increasing this. The two stage method is not suitable for this. \n     'None Shall Pass!' - The Black Knight")
    } else
      stop("first.stage.threshold is 0, consider increasing this. The two stage method is not suitable for this. \n     'You Shall Not Pass!' - Gandalf the Grey")
  }

  if(multicore != FALSE) stop("Package currently does not support parallel processing")

  if(abs(report.lowest.amount - round(report.lowest.amount)) > .Machine$double.eps^0.5 || report.lowest.amount < 1)
    warning("report.lowest.amount must be non-negative integer. Rounding up to non-negative integer")


  if(progress == 0) progress = FALSE

  this.call <- match.call()

  start.time <- proc.time()[3]
  if(multicore == FALSE){
    ts.output <- singlecore.twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold, progress)
  }

  if(!return.raw){
    ts.output$second.stage.sparse.matrix@x <- stats::p.adjust(ts.output$second.stage.sparse.matrix@x, method = multiple.hypotheses.correction)
  }


  report.lowest.amount = min(report.lowest.amount, length(ts.output$second.stage.sparse.matrix@x))

  fifth.lowest <- stats::quantile(ts.output$second.stage.sparse.matrix@x, probs = report.lowest.amount/length(ts.output$second.stage.sparse.matrix@x))
  indices.lowest.five <- Matrix::which(ts.output$second.stage.sparse.matrix <= fifth.lowest &
                                         ts.output$second.stage.sparse.matrix > 0 &
                                         ts.output$second.stage.sparse.matrix < 1, arr.ind = TRUE)
  lowest.five <- ts.output$second.stage.sparse.matrix[indices.lowest.five]
  attr(lowest.five, "names") = c(apply(Matrix::which(ts.output$second.stage.sparse.matrix < fifth.lowest &
                                                       ts.output$second.stage.sparse.matrix > 0, arr.ind = TRUE),
                                      1, function(.) paste0(., collapse = " x ")))

  total.runtime <- proc.time()[3] - start.time

  return.object <- list(most.significant.results = sort(lowest.five),
                        p.value.matrix = ts.output$second.stage.sparse.matrix,
                        marginal.significant = ts.output$passed.indices,
                        first.stage = ts.output$first.stage.p.values,
                        fst = first.stage.threshold,
                        runtime = total.runtime,
                        call = this.call)
  class(return.object) <- "twostageGWAS"

  return.object

}


#' Two stage method using single core
#'
#' @param survival.dataset the dataset
#' @param covariate.matrix the covariates
#' @param first.stage.threshold the fst
#'
#' @return a sparseMatrix object
#'
#'
#' @examples
#' print("make exmple")
singlecore.twostagecoxph <- function(survival.dataset, covariate.matrix, first.stage.threshold,
                                     progress = 50){
  first.stage.result <- firststagecoxph(survival.dataset, covariate.matrix, progress)

  first.stage.rejections <- (first.stage.result < first.stage.threshold)
  passed.indices <- which(first.stage.rejections == TRUE)
  amount.rejections <- max(sum(first.stage.rejections, na.rm = TRUE), 0, na.rm = TRUE)

  if(amount.rejections == 0) stop("No rejections in the first stage")
  if(amount.rejections == 1) stop("Only one rejection in the first stage")

  relevant.indices <- utils::combn(passed.indices, 2)
  second.stage.sparse.matrix <- Matrix::sparseMatrix(i = relevant.indices[1,], relevant.indices[2,],
                                             x = 1, triangular = TRUE)
  for (first.index in 1:(amount.rejections-1)){
    for (second.index in (first.index+1):amount.rejections){
      index.first.covariate  <- passed.indices[first.index]
      index.second.covariate <- passed.indices[second.index]
      covariate.one <- covariate.matrix[, index.first.covariate]
      covariate.two <- covariate.matrix[, index.second.covariate]

      fitted.model <- survival::coxph(survival.dataset ~ covariate.one * covariate.two)
      second.stage.sparse.matrix[index.first.covariate, index.second.covariate] <-
        summary(fitted.model)$coefficients[3,5]

      if(max((first.index*(amount.rejections - 1) + second.index) %% progress == 0, FALSE, na.rm = TRUE)){
        progress.frac <- (first.index*(amount.rejections - 1) + second.index)/(amount.rejections*(amount.rejections-1)/2)
        cat("\rSecond stage is at ", round(progress.frac*100, digits = 0), "% progress. ",
            "Estimated time until completion: ",
            round((1-progress.frac)/progress.frac*(proc.time()[3] - start.time.first.stage)/3600, digits = 2),
            " hours. \t\t\t\t", sep = "")
        utils::flush.console()
      }
    }
  }

  return(list(second.stage.sparse.matrix = second.stage.sparse.matrix, passed.indices = passed.indices,
              first.stage.p.values = first.stage.result))
}



firststagecoxph <- function(survival.dataset, covariate.matrix, progress = 50){
  start.time.first.stage <- proc.time()[3]
  p.value.vector <- rep(1, length = dim(covariate.matrix)[2])
  for(covariate.index in 1:length(p.value.vector)){
    this.covariate <- covariate.matrix[, covariate.index]
    fitted.model <- survival::coxph(survival.dataset ~ this.covariate)
    p.value.vector[covariate.index] <- summary(fitted.model)$coefficients[1,5]

    if(max(covariate.index %% progress == 0, FALSE, na.rm = TRUE)){
      progress.frac <- covariate.index/length(p.value.vector)
      cat("\r", "First stage is at ", round(progress.frac*100, digits = 0), "% progress. ",
          "Estimated time until completion: ",
          round((1-progress.frac)/progress.frac*(proc.time()[3] - start.time.first.stage)/60, digits = 1),
          " minutes. \t\t\t", sep = "")
      utils::flush.console()
    }
  }
  if(progress != 0) cat("\rFirst stage complete. Commencing second stage.")
  return(p.value.vector)
}

#' Print method for twostage object
#'
#' @param twostageGWASobject object of class twostageGWAS
#' @param ... optional arguments to print methods
#'
#' @return The twostageGWAS object invisibly
#' @export
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
#' twostagecoxph(survival.dataset, covariate.matrix, progress = 0)
print.twostageGWAS <- function(twostageGWASobject, ...){
  cat("\nCall:\n",
    paste(deparse(twostageGWASobject$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat(length(twostageGWASobject$marginal.significant), " covariates were marginally significant at level ",
      twostageGWASobject$fst, ", resulting in ", round(length(twostageGWASobject$marginal.significant)*(length(twostageGWASobject$marginal.significant)-1)/2),
      " second stage tests.\n\n", sep = "")
  if(length(twostageGWASobject$most.significant.results) == 0){
    cat("No significant results were found.")
  } else{
    cat("These are the most significant interactions found:\n")
        print(round(twostageGWASobject$most.significant.results[1:min(length(twostageGWASobject$most.significant.results), 5)],
                    digits = min(6L, getOption("digits"))))
  }
  invisible(twostageGWASobject)
}
