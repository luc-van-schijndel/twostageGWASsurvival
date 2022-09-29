#' (WIP) Perform a two stage analysis on a survival dataset to detect epistasis
#'
#' Work in progress. Only basic functionality works.
#'
#' @param survival.dataset The survival dataset describing the outcome.
#' @param covariate.matrix The nxp-matrix of covariates of the p covariates of the n patients.
#' @param first.stage.threshold numeric scalar denoting the threshold for the first stage. If a covariate
#'          has a p-value lower than this threshold, it will be passed on to the second stage.
#' @param multiple.hypotheses.correction Correction method, a character string. Passed to stats::p.adjust
#' @param multicore logical, default FALSE; whether or not the function should use multiple cores
#'                    in its calculations. See Details.
#' @param updatefile path to a text file where updates may be written. Necessary for parallel
#'                     computations, since the connection to the terminal will be lost. This
#'                     file will serve as a stand in for the terminal.
#' @param control Object of class twostagecoxph.control specifying various options for performance
#'                  of the two stage method. Default is twostagecoxph.control()
#'
#' @return A twostageGWAS object which is a list of 7 entries:
#'   \item{most.significant.results}{A list describing the most significant results found. The
#'   number of results reported is specified by the control parameter, default 5. The
#'   list contains 3 items:
#'   \describe{
#'     \item{\code{interacting.snps}}{the names of the interaction, in "name_snp_1 x name_snp_2" format.
#'     In case the covariates in \code{covariate.matrix} have a names attribute, these names are used.
#'     Otherwise the index of the covariates within the matrix is used.}
#'     \item{\code{p.value.epistasis}}{the corresponding p-values of the interactions. Unless \code{return.raw = TRUE}
#'     is specified in the control parameter, these p-values will be corrected for the multiple hypotheses
#'     tested with the method specified by the \code{multiple.hypotheses.correction} parameter.}
#'     \item{\code{duplicate.interactions}}{the list of duplicate interactions found, corresponding to
#'     the interactions specified in this list. An interaction is said to be duplicate if the corresponding
#'     p-value is the same up until 7 significant figures. }
#'   }}
#'   \item{p.value.matrix}{A \code{sparseMatrix} from the package "Matrix", specifying the resulting
#'   upper triangular p-value matrix obtained from the second stage.  Unless \code{return.raw = TRUE}
#'     is specified in the control parameter, these p-values will be corrected for the multiple hypotheses
#'     tested with the method specified by the \code{multiple.hypotheses.correction} parameter.}
#'   \item{marginal.significant}{A vector of named integers, specifying the indices of covariates which
#'   were found to be marginally significant in the first stage. }
#'   \item{first.stage}{A vector specifying the p-values found in the first stage. These p-values
#'   are \emph{not} corrected for the multiple tested hypotheses. }
#'   \item{fst}{The threshold for significance used in the first stage as specified in \code{first.stage.threshold}}
#'   \item{runtime}{The total runtime of the function in seconds.}
#'   \item{call}{The matched call.}
#' @export
#'
#' @note Parallel processing requires a properly registered parallel back-end, such as one obtained
#'         from "doParallel::registerDoParallel(2)" to be used by \code{foreach} and the \code{%dopar%} binary operator.
#'
#' @seealso \code{\link{foreach}}, \code{\link{print.twostageGWAS}}, \code{\link{twostagecoxph.control}}
#'
#' @importFrom survival coxph
#' @importFrom foreach %dopar%
#' @importFrom stats cor
#' @importFrom stats var
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
#' twostagecoxph(survival.dataset, covariate.matrix, control = twostagecoxph.control(progress = 0))
#'
#' str(example_survival_data)
#' str(example_snp_data)
#'
#' print(foo <- twostagecoxph(example_survival_data, example_snp_data[,1:300],
#'                             first.stage.threshold = 1e-5))
#' print(bar <- twostagecoxph(example_survival_data, example_snp_data[,1:300],
#'                             first.stage.threshold = 1e-4))
#' #[,1:300] subsetting is added to speed up the example. Try removing it! :)
#'
#' # As we can see, foo and bar have different results. A lower FST generally gives more power, but it
#' # it risks the possibility to be too strict and consequently *decreasing* power.
twostagecoxph <- function(survival.dataset, covariate.matrix, first.stage.threshold = 0.05,
                          multiple.hypotheses.correction = "bonferroni", multicore = FALSE,
                          updatefile = "", control = twostagecoxph.control()){
  report.lowest.amount = control$report.lowest.amount
  return.raw = control$return.raw
  progress = control$progress
  max.coef = control$max.coef
  max.batchsize = control$max.batchsize
  upper.bound.correlation = control$upper.bound.correlation

  # first.stage.threshold = 0.05; multiple.hypotheses.correction = "bonferroni"; multicore = FALSE; report.lowest.amount = 5; return.raw = FALSE; progress = 0; max.coef = 5; max.batchsize = 1000; updatefile = ""; upper.bound.correlation = 0.95
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

  if(max.batchsize < 2) stop("max.batchsize must be 2 or greater.")


  if(progress == 0) progress = FALSE

  this.call <- match.call()

  snps.are.named = !is.null(dimnames(covariate.matrix)[[2]])

  start.time <- proc.time()[3]
  if(multicore == FALSE){
    ts.output <- singlecore.twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold,
                                          progress, max.coef, upper.bound.correlation)
  }

  if(multicore != FALSE){
    ts.output <- multicore.twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold,
                                         progress, max.coef, max.batchsize, updatefile, upper.bound.correlation)
  }

  if(length(ts.output$passed.indices) < 2){
    total.runtime <- proc.time()[3] - start.time
    names(total.runtime) = c("seconds")

    return.object <- list(most.significant.results = NA,
                          p.value.matrix = NA,
                          marginal.significant = ts.output$passed.indices,
                          first.stage = ts.output$first.stage.p.values,
                          fst = first.stage.threshold,
                          runtime = total.runtime,
                          call = this.call)
    class(return.object) <- "twostageGWAS"



    return(return.object)
  }

  if(!return.raw){
    ts.output$second.stage.sparse.matrix@x <- stats::p.adjust(ts.output$second.stage.sparse.matrix@x, method = multiple.hypotheses.correction)
  }

  unique.p.values <- unique(ts.output$second.stage.sparse.matrix@x)

  report.lowest.amount = min(report.lowest.amount, sum(unique.p.values < 1, na.rm = TRUE))
  if(report.lowest.amount > 0){
    fifth.lowest.unique <- stats::quantile(unique.p.values,
                                    probs = report.lowest.amount/sum(!is.na(unique.p.values)),
                                    na.rm = TRUE)
    indices.lowest.five <- Matrix::which(ts.output$second.stage.sparse.matrix <= fifth.lowest.unique &
                                           ts.output$second.stage.sparse.matrix > 0 &
                                           ts.output$second.stage.sparse.matrix < 1, arr.ind = TRUE)

    lowest.five.dupl <- ts.output$second.stage.sparse.matrix[indices.lowest.five]
    lowest.table <- table(lowest.five.dupl) - 1
    if(snps.are.named){
      names.lowest.five <- matrix(c(dimnames(covariate.matrix)[[2]][as.vector(indices.lowest.five)]),
                                  ncol = 2)
      names(lowest.five.dupl) <- c(apply(names.lowest.five,
                                         1, function(.) paste0(., collapse = " x ")))
    } else{
      names(lowest.five.dupl) <- c(apply(indices.lowest.five,
                                         1, function(.) paste0(., collapse = " x ")))
    }


    lowest.five <- lowest.five.dupl[!duplicated(lowest.five.dupl)] #this keeps the "names" attribute

    lowest.five <- sort(lowest.five, decreasing = FALSE)
    attr(lowest.five, "duplicate results") = as.integer(table(lowest.five.dupl) - 1) #the table was already sorted



    duplicate.list <- vector("list", report.lowest.amount)
    names(duplicate.list) <- names(lowest.five)
    for(interaction in 1:report.lowest.amount){
      this.list = list(lowest.five.dupl[which(lowest.five.dupl == lowest.five[interaction])])
      names(this.list) = names(lowest.five)[interaction]
      duplicate.list[[interaction]] <- this.list[[1]][names(this.list[[1]]) != names(this.list)]
    }

    lowest.five.list <- list(interacting.snps = attr(lowest.five, "names"),
                             p.value.epistasis = as.numeric(sort(unique(lowest.five.dupl))),
                             duplicate.interactions = duplicate.list)
  } else lowest.five.list = list()
  total.runtime <- proc.time()[3] - start.time
  names(total.runtime) = c("seconds")

  clear.current.line()
  cat("\rAnalysis completed. Runtime:", total.runtime)

  return.object <- list(most.significant.results = lowest.five.list,
                        p.value.matrix = ts.output$second.stage.sparse.matrix,
                        marginal.significant = ts.output$passed.indices,
                        first.stage = ts.output$first.stage.p.values,
                        fst = first.stage.threshold,
                        runtime = total.runtime,
                        call = this.call)
  class(return.object) <- "twostageGWAS"

  return.object

}

#' Ancillary arguments for controlling the operation of the two stage method using twostagecoxph
#'
#' @description This is used to set various numeric parameters controlling the operation of the
#'                two stage method performed by \code{twostagecoxph}. Typically it would only be used
#'                in a call to \code{twostagecoxph}.
#'
#' @param report.lowest.amount integer, default 5; denotes how many of the most significant interactions
#'                               the function should report. Altering this value does
#'                               not affect the maximum reported by the print.twostageGWAS function
#' @param return.raw logical, default FALSE; whether or not the output should contain the raw p-values or
#'                     the multiple hypotheses corrected p-values.
#' @param progress numeric, default 1000;  how many iterations should pass silently until an update is given
#'                   about runtime and progress until completion of stages. Set to 0 for no output.
#' @param max.coef numeric, default 5; maximum value for all coefficients in the fitted models. If any
#'                   are larger than this (in absolute value), then the model rejected and ignored in
#'                   further analysis.
#' @param upper.bound.correlation numeric; the upper bound on the correlation between two
#'                                  covariates. If exceeded (in absolute value), the model
#'                                  is not fitted.
#' @param max.batchsize maximum size of one batch, default 1000; for parallel computing,
#'                        should be lowered if issues arise concerning memory. Can be raised
#'                        to slightly increase performance
#'
#' @return a list containing the values of each of the above constants.
#' @export
#'
twostagecoxph.control <- function(report.lowest.amount = 5, return.raw = FALSE, progress = 1000,
                                  max.coef = 5, max.batchsize = 1000, upper.bound.correlation = 0.9){
  if((!is.numeric(report.lowest.amount)) ||
     (trunc(report.lowest.amount) != report.lowest.amount)[1] ||
     (length(report.lowest.amount) != 1)) stop("report.lowest.amount must be an integer scalar")
  if(!is.numeric(progress) ||
     (trunc(progress) != progress)[1] ||
     length(progress) != 1) stop("progress must be an integer scalar")
  if(!is.numeric(max.coef) ||
     (trunc(max.coef) != max.coef)[1] ||
     length(max.coef) != 1) stop("max.coef must be an integer scalar")
  if(!is.numeric(max.batchsize) ||
     (trunc(max.batchsize) != max.batchsize)[1] ||
     length(max.batchsize) != 1) stop("max.batchsize must be an integer scalar")
  if(!is.numeric(upper.bound.correlation) ||
     (upper.bound.correlation[1] < 0 || upper.bound.correlation[1] > 1) ||
     length(upper.bound.correlation) != 1) stop("upper.bound.correlation must be an double scalar in the interval [0,1]")
  if(!is.logical(return.raw) || length(return.raw) != 1) stop("return.raw must be a logical scalar")
  if(max.batchsize < 2) stop("max.batchsize must be 2 or greater")

  return(list(report.lowest.amount = report.lowest.amount, return.raw = return.raw, progress = progress,
              max.coef = max.coef, max.batchsize = max.batchsize, upper.bound.correlation = upper.bound.correlation))
}

prefitting.check.one <- function(covariate) {
  passed <- TRUE
  if (stats::var(covariate, na.rm = TRUE) < 0.1)
    passed <- FALSE
  return(passed)
}

prefitting.check.two <- function(covariate.one, covariate.two, upper.bound.correlation) {
  passed <- TRUE
  if (max(stats::cor(matrix(
    c(covariate.one, covariate.two, covariate.one * covariate.two),
    ncol = 3
  ))[upper.tri(matrix(0, 3, 3))],
  na.rm = TRUE) > upper.bound.correlation)
    passed <- FALSE
  return(passed)
}


convergence.check <- function(coxph.model, max.coef = 5){
  converged <- TRUE
  if(max(abs(coxph.model$coefficients), na.rm = TRUE) > max.coef) converged <- FALSE
  if(coxph.model$iter >= 20) FALSE
  return(converged)
}


#' Calculates the optimal batch configuration
#'
#' @param no.covariates number of covariates to distribute
#' @param max.batchsize upper limit on amount covariates in each batch
#' @param no.workers number of workers
#'
#' @return a list of size 3, containing the number of covariates in each batch, the number of batches
#' (including possible empty ones at the end, this way it always is a multiple of no.workers),
#' and a vector of batchsizes for the last iteration (including trailing zeros)
#' @export
#'
#' @details configures this in such a way that the number of covariates does not exceed the
#' allowed maximum, while having each batch filled equally (not necessarily full). The last batch
#' may be not fully filled. Made for the first stage, but easily adaptible to the second stage.
optimal.batch.configuration <- function(no.covariates, max.batchsize, no.workers = 1){
  intermediate.no.batches <- min(ceiling(no.covariates/no.workers), max.batchsize)
  intermediate.batchsize <- ceiling(no.covariates/intermediate.no.batches)
  final.iters <- ceiling(intermediate.batchsize/no.workers)
  final.no.batches <- final.iters*no.workers
  final.batchsize <- ceiling(no.covariates/final.no.batches)

  last.iter.covs <- no.covariates - final.batchsize*(final.iters-1)*no.workers
  last.batchsize <- ceiling(last.iter.covs/no.workers)
  last.no.full.batches <- ceiling(last.iter.covs/last.batchsize - 1)
  last.partial.batch <- last.iter.covs - last.batchsize*last.no.full.batches
  last.batchsizes <- c(rep(last.batchsize, last.no.full.batches), last.partial.batch)
  #add trailing 0's
  last.batchsizes <- c(last.batchsizes, rep(0, no.workers - length(last.batchsizes)))

  return(list(optimal.batchsize = final.batchsize, optimal.no.batches = final.no.batches,
              last.batchsizes = last.batchsizes))
}

singlecore.twostagecoxph <- function(survival.dataset, covariate.matrix, first.stage.threshold,
                                     progress = 50, max.coef = 5, upper.bound.correlation = 0.95){
  first.stage.result <- firststagecoxph(survival.dataset, covariate.matrix, progress, max.coef)

  first.stage.rejections <- (first.stage.result < first.stage.threshold)
  passed.indices <- which(first.stage.rejections == TRUE)
  amount.rejections <- max(sum(first.stage.rejections, na.rm = TRUE), 0, na.rm = TRUE)

  if(amount.rejections == 0) warning("No rejections in the first stage")
  if(amount.rejections == 1) warning("Only one rejection in the first stage")
  if(amount.rejections < 2) return(list(second.stage.sparse.matrix = NA, passed.indices = passed.indices,
                                        first.stage.p.values = first.stage.result))

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
      if(prefitting.check.two(covariate.one, covariate.two, upper.bound.correlation)) {
        tryCatch(
          fitted.model <-
            coxph(survival.dataset ~ covariate.one * covariate.two),
          warning = function(w) {
            if (grepl("coefficient may be infinite", w$message)) {
              #print("An error was given, which is taken into account in convergence.check")
            } else if (grepl("out of iterations", w$message)) {
              #print("An error was given, which is taken into account in convergence.check")
            }
            else {
              message(w$message)
            }
          }
        )



        if (convergence.check(fitted.model, max.coef)) {
          second.stage.sparse.matrix[index.first.covariate, index.second.covariate] <-
            summary(fitted.model)$coefficients[3, 5]
        }
      } else {
        second.stage.sparse.matrix[index.first.covariate, index.second.covariate] <-
          NA
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



#' Multicore method of performing the second stage
#'
#' @param survival.dataset the outcome data
#' @param covariate.matrix The SNPs
#' @param first.stage.threshold the FST
#' @param progress set to 0 for no updates
#' @param max.coef maximum value of fitted weights before declared non-converged
#' @param updatefile path to that file
#' @param max.batchsize max number of covariates in one batch
#' @param upper.bound.correlation upper bound on the correlation before not checked
#'
#' @return list of p-value matrix, first stage p-values and which ones passed.
#'
#' @details Similarly to the multicore method of the first stage, this function works with
#'   batches of covariates to alleviate possible memory issues. The optimal size of the batches
#'   is calculated in a similar fashion as during the first stage, only here we halve the maximum
#'   batchsize, since (almost always) two batches of covariates will be in memory at the same time.
#'
#'   The testing for interactions is done in a first-in, last-out approach. The first batch of
#'   covariates will be tested for interactions with itself, then for with all covariates from
#'   subsequent batches. The second batch does not need to test for interactions with the first one,
#'   since the first one already did that. This allows the second batch to be done before the first one,
#'   hence the "first-in, last-out" naming. The number of batches will always be a multiple of
#'   two times the number of worker cores; this should ensure that all workers should be done
#'   at the same time.
#'
#'   There is quite some code-duplication going on (DRY! OMG! U NUB), but that's for a reason. Function
#'   calls in R have a tiny bit of overhead, and since we would do these function-calls a lot (and I mean a LOT)
#'   of times, I have opted to simply copy the function's code into this function itself, and forego
#'   the function itself.
#'
#'
#'   Since this function is not exported, I will discusse the 6 indices found in the for-loops
#'   (at some point we are 4 for-loops deep)
#'   first.batch.real.indices
#'   first/second.local.index
#'   first/second.covariate.index
#'
#'   passed.indices are not used, since amount rejections is better since possible NAs have been
#'   removed.
#'
multicore.twostagecoxph <- function(survival.dataset, covariate.matrix, first.stage.threshold,
                                     progress = 50, max.coef = 5, updatefile = "", max.batchsize = 1000,
                                    upper.bound.correlation = 0.95){
  first.stage.result <- firststagecoxph.multicore(survival.dataset, covariate.matrix, progress,
                                                  max.coef, max.batchsize = 1000, updatefile)

  #trailing NAs have names "", so we drop these
  first.stage.result <- first.stage.result[names(first.stage.result) != ""]

  #revert (possible) reordering introduced in first stage
  first.stage.permutation <- order(as.integer(names(first.stage.result)))
  first.stage.result <- first.stage.result[first.stage.permutation]

  first.stage.rejections <- (first.stage.result < first.stage.threshold)
  passed.indices <- which(first.stage.rejections == TRUE)
  amount.rejections <- max(sum(first.stage.rejections, na.rm = TRUE), 0, na.rm = TRUE)

  if(amount.rejections == 0) warning("No rejections in the first stage")
  if(amount.rejections == 1) warning("Only one rejection in the first stage")
  if(amount.rejections < 2) return(list(second.stage.sparse.matrix = NA, passed.indices = passed.indices,
                                        first.stage.p.values = first.stage.result))

  relevant.indices <- utils::combn(passed.indices, 2)
  second.stage.sparse.matrix <- Matrix::sparseMatrix(i = relevant.indices[1,], relevant.indices[2,],
                                                     x = 1, triangular = TRUE)

  no.workers <- foreach::getDoParWorkers()
  optimal.batchconf <- optimal.batch.configuration(no.covariates = ceiling(amount.rejections/2),
                                                   max.batchsize = floor(max.batchsize/2),
                                                   no.workers = no.workers)
  optimal.batchsize <- optimal.batchconf$optimal.batchsize
  optimal.no.batches <- optimal.batchconf$optimal.no.batches
  #$last.batchsizes

  #round up in the first half, down in the second
  batch.indices.matrix.first.half <- matrix(c(seq_len(sum(ceiling(amount.rejections/2))),
                                              rep(NA, optimal.batchsize*optimal.no.batches -
                                                    sum(ceiling(amount.rejections/2)))),
                                            nrow = optimal.batchsize, ncol = optimal.no.batches)
  batch.indices.matrix.second.half <- matrix(c(seq_len(sum(floor(amount.rejections/2))),
                                              rep(NA, optimal.batchsize*optimal.no.batches -
                                                    sum(floor(amount.rejections/2)))),
                                            nrow = optimal.batchsize, ncol = optimal.no.batches) +
                                      ceiling(amount.rejections/2)
  batch.indices.matrix <- cbind(batch.indices.matrix.first.half, batch.indices.matrix.second.half)
  rm(batch.indices.matrix.first.half, batch.indices.matrix.second.half)

  start.time.second.stage <- proc.time()[3]

  first.batch.index <- NULL #to suppress a note from devtools::check()

  output.matrix <- foreach::foreach(first.batch.index = 1:(dim(batch.indices.matrix)[2]), #no '-1' here, since the last batch must be checked with itself, contrary to covariates where this is not the case
                                    .packages = c("survival", "utils"),
                                    .export = c("prefitting.check.two", "convergence.check",
                                                "clear.current.line"),
                                    .combine = rbind) %dopar% {
    #I'm NOT gonna reindent these lines to match upt with the foreach-loop...
    first.batch.real.indices <- batch.indices.matrix[,first.batch.index]
    first.batch.real.indices <- first.batch.real.indices[!is.na(first.batch.real.indices)]

    #prepare the matrix we will return (for this first.batch)
    return.matrix <- matrix(NA, nrow = length(first.batch.real.indices), ncol = amount.rejections)

    #edge case: if this batch is the last one, and all indices are empty (can occur if batchsize = 1)
    if(length(first.batch.real.indices) == 0) return(return.matrix)

    #first check the batch with itself, if there is more than 1 covariate in the batch
    if(length(first.batch.real.indices) > 1){
    for(first.local.index in 1:(length(first.batch.real.indices)-1)){
      for(second.local.index in (first.local.index + 1):length(first.batch.real.indices)){
        index.first.covariate  <- first.batch.real.indices[first.local.index]
        index.second.covariate <- first.batch.real.indices[second.local.index]
        covariate.one <- covariate.matrix[, index.first.covariate]
        covariate.two <- covariate.matrix[, index.second.covariate]
        if(prefitting.check.two(covariate.one, covariate.two, upper.bound.correlation)){
          tryCatch(
            fitted.model <-
              coxph(survival.dataset ~ covariate.one * covariate.two),
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
        }


        if(convergence.check(fitted.model, max.coef)){
          return.matrix[first.local.index, second.local.index] <-
            summary(fitted.model)$coefficients[3,5]
        }

        if(max((amount.rejections*(first.local.index-1) + second.local.index - first.local.index*(first.local.index+1)/2) %% progress == 0, FALSE, na.rm = TRUE)){
          progress.frac <- (amount.rejections*(first.local.index-1) + second.local.index - first.local.index*(first.local.index+1)/2)/(amount.rejections*(amount.rejections-1)/2)

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
    }

    #now we check the first batch with (possible) second batches, if this is not the last batch
    if(dim(batch.indices.matrix)[2] != first.batch.index){
      for(second.batch.index in (first.batch.index + 1):dim(batch.indices.matrix)[2]){
        second.batch.real.indices <- batch.indices.matrix[,second.batch.index]
        second.batch.real.indices <- second.batch.real.indices[!is.na(second.batch.real.indices)]

        if(length(second.batch.real.indices) > 0){ #only if the second batch is not empty
        for(first.local.index in 1:(length(first.batch.real.indices))){
          for(second.local.index in 1:length(second.batch.real.indices)){
            index.first.covariate  <- first.batch.real.indices[first.local.index]
            index.second.covariate <- second.batch.real.indices[second.local.index]
            covariate.one <- covariate.matrix[, index.first.covariate]
            covariate.two <- covariate.matrix[, index.second.covariate]
            if(prefitting.check.two(covariate.one, covariate.two, upper.bound.correlation)){
              tryCatch(
                fitted.model <-
                  coxph(survival.dataset ~ covariate.one * covariate.two),
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
            }


            if(convergence.check(fitted.model, max.coef)){
              return.matrix[first.local.index, index.second.covariate] <- #note the discrepancy between the two subset-indices. This is intentional!
                summary(fitted.model)$coefficients[3,5]
            }

            if(max((amount.rejections*(first.local.index-1) + second.local.index - first.local.index*(first.local.index+1)/2) %% progress == 0, FALSE, na.rm = TRUE)){
              progress.frac <- (amount.rejections*(first.local.index-1) + second.local.index - first.local.index*(first.local.index+1)/2)/(amount.rejections*(amount.rejections-1)/2)

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
        }
      }
    }
    return(return.matrix)
  }

  #output.matrix is a dense matrix, which we need to expand to a sparse matrix with proper indexing
  non.na.indices <- which(!is.na(output.matrix), arr.ind = TRUE)
  second.stage.sparse.matrix <- Matrix::sparseMatrix(i = non.na.indices[,1], j = non.na.indices[,2],
                                                     x = as.vector(output.matrix[non.na.indices]), triangular = TRUE)

  return(list(second.stage.sparse.matrix = second.stage.sparse.matrix, passed.indices = passed.indices,
              first.stage.p.values = first.stage.result))
}


firststagecoxph <- function(survival.dataset, covariate.matrix, progress = 50, max.coef = 5){
  start.time.first.stage <- proc.time()[3]
  p.value.vector <- rep(1, length = dim(covariate.matrix)[2])
  for(covariate.index in 1:length(p.value.vector)){
    this.covariate <- covariate.matrix[, covariate.index]
    if(prefitting.check.one(this.covariate)){
      tryCatch(
        fitted.model <- coxph(survival.dataset ~ this.covariate),
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
    }
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
  names(p.value.vector) <- 1:length(p.value.vector)
  return(p.value.vector)
}


#' Perform the first stage of a GWAS in parallel
#'
#' @param survival.dataset the outcomes
#' @param covariate.matrix the covariates
#' @param progress how often progress must be reported
#' @param max.coef what the maximum coefficients may be when fitting
#' @param max.batchsize what the maximum batchsize must be.
#' @param updatefile path to file to replace terminal
#'
#' @return the p-values of the first stage (in random order) with proper names() attribute.
#'          As a consequence of possible empty batches, there may be trailing NA's. Ironically however,
#'          these are  not guaranteed to be trailing, due to the fact that batches can finish
#'          in any order. They are identified by their names attribute being empty: "".
#'
#' @importFrom foreach %dopar%
#' @importFrom survival coxph
#'
#' @details
#' The covariates are split into a number of batches. The amount of covariates in every batch
#' is optimal, i.e. the load is distributed equally over all batches, the number
#' of covariates in each batch does not exceed max.batchsize, and the number of batches is a
#' multiple of the number of worker-processes.
#' These batches are then analysed in parallel by the worker-processes that return the
#' corresponding p-values of the marginal associations of the covariates. These p-values have
#' their names attribute set to the index of the covariate, so that the processes do not need
#' to finish in order. If they do not, the order of the returned p-values will not match up
#' with the order in which the covariates were provided, so this problem is avoided.
#'
#'
#'
firststagecoxph.multicore <- function(survival.dataset, covariate.matrix, progress = 50,
                                      max.coef = 5, max.batchsize = 1000, updatefile = ""){
  start.time.first.stage <- proc.time()[3]

  no.covariates <- dim(covariate.matrix)[2]
  no.workers <- foreach::getDoParWorkers()
  optimal.batchconf <- optimal.batch.configuration(no.covariates, max.batchsize, no.workers)

  optimal.batchsize <- optimal.batchconf$optimal.batchsize
  no.processes <- optimal.batchconf$optimal.no.batches
  batchsizes.vector <- c(rep(optimal.batchsize, no.processes - no.workers),
                         optimal.batchconf$last.batchsizes)
  matrix.indices.processes <- matrix(c(seq_len(sum(batchsizes.vector)),
                                       rep(NA, no.workers*max(optimal.batchconf$last.batchsizes) -
                                             sum(optimal.batchconf$last.batchsizes))),
                                     ncol = length(batchsizes.vector),
                                     nrow = optimal.batchsize)

  process.index <- NULL #suppressing a note from devtools::check()
  p.value.vector <- foreach::foreach(process.index = seq_len(no.processes),
                                     .packages = c("survival", "utils"),
                                     .export = c("prefitting.check.one", "convergence.check",
                                                 "clear.current.line"),
                                     .combine = c) %dopar% {
    this.process.indices = (matrix.indices.processes[,process.index])[!is.na(matrix.indices.processes[,process.index])]
    if(length(this.process.indices) == 0) return(rep(NA, optimal.batchsize)) #if the batch is empty, return NA's
    this.process.p.values <- rep(NA, length(this.process.indices))
    for(covariate.index in 1:length(this.process.indices)){
      this.covariate <- covariate.matrix[, this.process.indices[covariate.index]]
      if(prefitting.check.one(this.covariate)){
        tryCatch(
          fitted.model <- survival::coxph(survival.dataset ~ this.covariate),
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
      }
      if (convergence.check(fitted.model, max.coef)) {
        this.process.p.values[covariate.index] <-
          summary(fitted.model)$coefficients[1, 5]
      }
      if(max(this.process.indices[covariate.index] %% progress == 0, FALSE, na.rm = TRUE)){
        progress.frac <- this.process.indices[covariate.index]/no.covariates
        clear.current.line()
        cat("\r", "First stage is at ", round(progress.frac*100, digits = 0), "% progress. ",
            "Estimated time until completion first stage: ",
            round((1-progress.frac)/progress.frac*(proc.time()[3] - start.time.first.stage)/60, digits = 1),
            " minutes. ", sep = "", file = updatefile)
        utils::flush.console()
      }
    }

    names(this.process.p.values) <- this.process.indices
    return(this.process.p.values)
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

#' Print method for twostageGWAS objects
#'
#' @usage \method{print}{twostageGWAS}(x, \ldots)
#' @param x object of class twostageGWAS
#' @param ... optional arguments passed on to cat and print.default functions.
#'
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
#' twostagecoxph(survival.dataset, covariate.matrix, control = twostagecoxph.control(progress = 0))
print.twostageGWAS <- function(x, ...){
  #if(class(x) != "twostageGWAS") stop("class must be twostageGWAS")
  clear.current.line()
  cat("\nCall:\n",
    paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "", ...)
  cat(length(x$marginal.significant), " covariates were marginally significant at level ",
      x$fst, ", resulting in ", round(length(x$marginal.significant)*(length(x$marginal.significant)-1)/2),
      " second stage tests.\n\n", sep = "", ...)
  if(length(x$most.significant.results) == 0){
    cat("No relevant results were found after applying the multiple hypotheses correction.", ...)
  } else{
    cat("These are the most significant interactions found with their respective p-values:\n", ...)
    #a maximum of 5 will be displayed
    similar.interactions <- rep(0, min(length(x$most.significant.results$duplicate.interactions), 5))
    for (i in 1:length(similar.interactions)){
      similar.interactions[i] <- length(x$most.significant.results$duplicate.interactions[[i]])
    }
    format.matrix <- matrix(c(x$most.significant.results$p.value.epistasis[1:length(similar.interactions)],
                              similar.interactions), nrow = 2, ncol = length(similar.interactions),
                            dimnames = list(c("p-value epistasis", "amount similar results"), x$most.significant.results$interacting.snps),
                            byrow = TRUE)
    print.default(format(format.matrix, drop0trailing = TRUE, digits = min(6L, getOption("digits"))), quote = FALSE, ...)
  }
  cat("\n", ...)
  invisible(x)
}
