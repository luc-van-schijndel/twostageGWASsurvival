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
#' @param max.coef numeric, default 5; maximum value for all coefficients in the fitted models. If any
#'                   are larger than this (in absolute value), then the model rejected and ignored in
#'                   further analysis.
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
#'
#' ## Not run:
#' #load("data/example_survival_data.rda")
#' #load("data/example_snp_data.rda")
#' #str(example_survival_data)
#' #str(example_snp_data)
#' #foo <- twostagecoxph(example_survival_data, example_snp_data, first.stage.threshold = 1e-5)
#' ## End(Not run)
twostagecoxph <- function(survival.dataset, covariate.matrix, first.stage.threshold = 0.05,
                          multiple.hypotheses.correction = "bonferroni", multicore = FALSE,
                          report.lowest.amount = 5, return.raw = FALSE,
                          progress = 1000, max.coef = 5){
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
    ts.output <- singlecore.twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold,
                                          progress, max.coef)
  }

  if(!return.raw){
    ts.output$second.stage.sparse.matrix@x <- stats::p.adjust(ts.output$second.stage.sparse.matrix@x, method = multiple.hypotheses.correction)
  }


  report.lowest.amount = min(report.lowest.amount, sum(ts.output$second.stage.sparse.matrix@x < 1))

  fifth.lowest <- stats::quantile(ts.output$second.stage.sparse.matrix@x,
                                  probs = report.lowest.amount/sum(!is.na(ts.output$second.stage.sparse.matrix@x)),
                                  na.rm = TRUE)
  indices.lowest.five <- Matrix::which(ts.output$second.stage.sparse.matrix <= fifth.lowest &
                                         ts.output$second.stage.sparse.matrix > 0 &
                                         ts.output$second.stage.sparse.matrix < 1, arr.ind = TRUE)
  lowest.five <- ts.output$second.stage.sparse.matrix[indices.lowest.five]
  attr(lowest.five, "names") = c(apply(Matrix::which(ts.output$second.stage.sparse.matrix <= fifth.lowest &
                                                       ts.output$second.stage.sparse.matrix > 0, arr.ind = TRUE),
                                      1, function(.) paste0(., collapse = " x ")))

  total.runtime <- proc.time()[3] - start.time
  names(total.runtime) = c("seconds")

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

convergence.check <- function(coxph.model, max.coef = 5){
  converged <- TRUE
  if(max(abs(coxph.model$coefficients), na.rm = TRUE) > max.coef) converged <- FALSE
  if(coxph.model$iter >= 20) FALSE
  return(converged)
}


singlecore.twostagecoxph <- function(survival.dataset, covariate.matrix, first.stage.threshold,
                                     progress = 50, max.coef = 5){
  first.stage.result <- firststagecoxph(survival.dataset, covariate.matrix, progress, max.coef)

  first.stage.rejections <- (first.stage.result < first.stage.threshold)
  passed.indices <- which(first.stage.rejections == TRUE)
  amount.rejections <- max(sum(first.stage.rejections, na.rm = TRUE), 0, na.rm = TRUE)

  if(amount.rejections == 0) stop("No rejections in the first stage")
  if(amount.rejections == 1) stop("Only one rejection in the first stage")

  relevant.indices <- utils::combn(passed.indices, 2)
  second.stage.sparse.matrix <- Matrix::sparseMatrix(i = relevant.indices[1,], relevant.indices[2,],
                                             x = 1, triangular = TRUE)
  start.time.second.stage <- proc.time()[3]
  for (first.index in 1:(amount.rejections-1)){
    for (second.index in (first.index+1):amount.rejections){
      index.first.covariate  <- passed.indices[first.index]
      index.second.covariate <- passed.indices[second.index]
      covariate.one <- covariate.matrix[, index.first.covariate]
      covariate.two <- covariate.matrix[, index.second.covariate]

      tryCatch(
        fitted.model <-
          survival::coxph(survival.dataset ~ covariate.one * covariate.two),
        warning = function(w) {
          if (grepl("coefficient may be infinite.", w$message)) {
            #print("An error was given, which is taken into account in convergence.check")
          } else if (grepl("out of iterations", w$message)) {
            #print("An error was given, which is taken into account in convergence.check")
          }
          else {
            message(w$message)
          }
        }
      )

      #fitted.model <- survival::coxph(survival.dataset ~ covariate.one * covariate.two)
      if(convergence.check(fitted.model, max.coef)){
        second.stage.sparse.matrix[index.first.covariate, index.second.covariate] <-
          summary(fitted.model)$coefficients[3,5]
      }

      if(max((amount.rejections*(first.index-1) + second.index - first.index*(first.index+1)/2) %% progress == 0, FALSE, na.rm = TRUE)){
        progress.frac <- (amount.rejections*(first.index-1) + second.index - first.index*(first.index+1)/2)/(amount.rejections*(amount.rejections-1)/2)

        estimated.seconds <- (1-progress.frac)/progress.frac*(proc.time()[3] - start.time.second.stage)
        text.estimated.time <- ""
        if(estimated.seconds >=3600) text.estimated.time = paste0(round(estimated.seconds/3600, digits = 1), " hours.")
        if(estimated.seconds < 3600) text.estimated.time = paste0(round(estimated.seconds/60, digits = 1), " minutes.")
        if(estimated.seconds < 60)   text.estimated.time = "less than a minute."
        clear.current.line()
        cat("\rSecond stage is at ", round(progress.frac*100, digits = 1), "% progress. ",
            "Estimated time until completion: ",
            text.estimated.time, sep = "")
        utils::flush.console()
      }
    }
  }

  return(list(second.stage.sparse.matrix = second.stage.sparse.matrix, passed.indices = passed.indices,
              first.stage.p.values = first.stage.result))
}



firststagecoxph <- function(survival.dataset, covariate.matrix, progress = 50, max.coef = 5){
  start.time.first.stage <- proc.time()[3]
  p.value.vector <- rep(1, length = dim(covariate.matrix)[2])
  for(covariate.index in 1:length(p.value.vector)){
    this.covariate <- covariate.matrix[, covariate.index]
    tryCatch(
      fitted.model <- survival::coxph(survival.dataset ~ this.covariate),
      warning = function(w) {
        if (grepl("coefficient may be infinite.", w$message)) {
          #print("An error was given, which is taken into account in convergence.check")
        } else {
          message(w$message)
        }
      }
    )
    fitted.model <-
      survival::coxph(survival.dataset ~ this.covariate)
    if (convergence.check(fitted.model, max.coef)) {
      p.value.vector[covariate.index] <-
        summary(fitted.model)$coefficients[1, 5]
    }
    if(max(covariate.index %% progress == 0, FALSE, na.rm = TRUE)){
      progress.frac <- covariate.index/length(p.value.vector)
      clear.current.line()
      cat("\r", "First stage is at ", round(progress.frac*100, digits = 0), "% progress. ",
          "Estimated time until completion first stage: ",
          round((1-progress.frac)/progress.frac*(proc.time()[3] - start.time.first.stage)/60, digits = 1),
          " minutes. ", sep = "")
      utils::flush.console()
    }
  }

  if(progress != 0){
    clear.current.line()
    cat("\rFirst stage complete. Commencing second stage. ")
  }
  return(p.value.vector)
}

clear.current.line <- function(){
  cat("\r", rep(" ", (getOption("width")+4)/2))
  utils::flush.console()
}

#' Print method for twostage object
#'
#' @param x object of class twostageGWAS
#' @param ... optional arguments passed on to cat and print.default functions.
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
print.twostageGWAS <- function(x, ...){
  cat("\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "", ...)
  cat(length(x$marginal.significant), " covariates were marginally significant at level ",
      x$fst, ", resulting in ", round(length(x$marginal.significant)*(length(x$marginal.significant)-1)/2),
      " second stage tests.\n\n", sep = "", ...)
  if(length(x$most.significant.results) == 0){
    cat("No relevant results were found after applying the multiple hypotheses correction.", ...)
  } else{
    cat("These are the most significant interactions found with their respective p-values:\n", ...)
        print.default(round(x$most.significant.results[1:min(length(x$most.significant.results), 5)],
                    digits = min(6L, getOption("digits"))), ...)
  }
  invisible(x)
}
