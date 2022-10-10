#' Performs a large scale two stage analysis on a survival dataset to detect epistasis using batches
#'
#' An adaptation of \code{\link{twostagecoxph}} where the covariates will be read from multiple files.
#' This is to alleviate possible memory constrains. The main functionality assumes a parallel back-end
#' is registered, e.g. using doParallel::registerDoParallel(2).
#'
#' @param survival.dataset The survival dataset describing the outcome.
#' @param covariate.filepaths The vector of paths of the files containing the covariates. See Details.
#' @param first.stage.threshold numeric scalar denoting the threshold for the first stage. If a covariate
#'          has a p-value lower than this threshold, it will be passed on to the second stage.
#' @param multiple.hypotheses.correction Correction method, a character string. Passed to \code{\link{stats::p.adjust}}
#' @param updatefile path to a text file where updates may be written. Necessary for parallel
#'                     computations, since the connection to the terminal will be lost. This
#'                     file will in that case serve as a stand-in for the terminal.
#' @param control object of class \code{twostagecoxph.control} specifying various options for performance
#'                  of the two stage method.
#' @param number.of.covariates the combined total of covariates found in the files specified.
#'                               If left unspecified, the function will automatically determine this
#'                               value by reading through all files before starting the analysis, thereby
#'                               increasing runtimes.
#' @param read.function function used to read batches of covariates from the files specified.
#'                        By default, this is \code{function(x) as.matrix(utils::read.table(x))}. Note that
#'                        the first argument of the function will be used to pass the filepaths. See Details.
#'
#' @details   At all times on all registered cores, the contents of at most 2 files will be in memory.
#'            Take note however, that the results will be written to a matrix of substantial size
#'            (\code{covariates_per_file x total_number_of_covariates}) so it is best if enough memory
#'            is still available after reading the files.
#'
#'            The function assumes that the object returned by \code{read.function} is a matrix containing
#'            the covariates with one column per covariate and one row per subject This assumption
#'            will only be checked once at function start. If the files are not structured the
#'            same, the results may be unreliable or an error may be thrown. If the \code{read.function}
#'            assignes dimnames to the columns, these names will be kept and passed on to the resulting
#'            object. If not, the indices of the covariates (counting across all files) will be used as
#'            the names.
#'
#'            The function
#'            assumes the files have same number of covariates, except for the last one,
#'            which may have less. \code{number.of.covariates}
#'            will be divided by the amount of covariates per file, and the remainder of that division
#'            \emph{must} be the number of covariates in the last file.
#'
#' @return A twostageGWAS object, which is a list of 7 entries:
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
#'     tested with the method specified by the \code{multiple.hypotheses.correction} parameter. The names
#'     of the dimensions of the matrix match the ones specified in the files, if the read.function assigns
#'     dimnames to the matrix. The row and columns corresponding to covariates that were not passed to
#'     the second stage, have names '"NA"' (Note: not 'NA'). The entries in these rows and columns are also all absent.}
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
#' @importFrom utils read.table
#'
#' @examples
#' str(example_survival_data)
#' str(example_snp_data)
#'
#' #Split the covariate matrix into various files.
#' number.of.covs <- dim(example_snp_data)[2]
#' number.of.files <- 6 #500 is not a multiple of 500, so the last file has less covariates than the other ones
#' temp.snpfile.paths <- tempfile(rep("snpfile", number.of.files),
#'                                tmpdir = tempdir(check = TRUE),
#'                                fileext = ".txt")
#' indices.matrix <- matrix(c(seq_len(number.of.covs),
#'                            rep(NA, ceiling(number.of.covs/number.of.files)*number.of.files - number.of.covs)),
#'                          ncol = number.of.files)
#' for(file.num in 1:number.of.files){
#'   indices <- indices.matrix[,file.num]
#'   write.table(example_snp_data[,indices[!is.na(indices)]], file = temp.snpfile.paths[file.num])
#' }
#'
#'
#' print(foo <- batched.twostagecoxph(example_survival_data, temp.snpfile.paths,
#'                                    number.of.covariates = number.of.covs,
#'                                    first.stage.threshold = 1e-2))
batched.twostagecoxph <- function(survival.dataset, covariate.filepaths, first.stage.threshold = 0.05,
                          multiple.hypotheses.correction = "bonferroni",
                          updatefile = "", control = twostagecoxph.control(),
                          number.of.covariates = 0,
                          read.function = function(x) as.matrix(read.table(x))){
  report.lowest.amount = control$report.lowest.amount
  return.raw = control$return.raw
  progress = control$progress
  max.coef = control$max.coef
  max.batchsize = control$max.batchsize
  upper.bound.correlation = control$upper.bound.correlation

  # first.stage.threshold = 0.05; multiple.hypotheses.correction = "bonferroni"; multicore = FALSE; report.lowest.amount = 5; return.raw = FALSE; progress = 0; max.coef = 5; max.batchsize = 84; updatefile = ""; upper.bound.correlation = 0.95; snps.are.named = FALSE
  if(!survival::is.Surv(survival.dataset)) stop("survival.dataset must be a Surv object")
  test.covariates <- read.function(covariate.filepaths[1])
  if(!is.matrix(test.covariates)) stop("covariates read from files are not properly structured")
  if(dim(test.covariates)[1] != dim(survival.dataset)[1]) stop("covariates read from files do not have same wrong dimensions")


  if(!(dim(survival.dataset)[1] %in% dim(test.covariates))) stop("No dimension of covariate.matrix matches dimension of survival.dataset")
  if(dim(survival.dataset)[1] != dim(test.covariates)[1] &&
     dim(survival.dataset)[1] == dim(test.covariates)[2]) {
    stop("The covariate matrices present in files are probably transposed")
  }


  if(!is.numeric(first.stage.threshold)) stop("first.stage.threshold must be numeric")
  if(length(first.stage.threshold) != 1) stop("first.stage.threshold must be of length 1")
  if(first.stage.threshold > 1 || first.stage.threshold < 0) stop("first.stage.threshold must be in the interval [0,1]")
  if(first.stage.threshold == 1) stop("first.stage.threshold is 1, consider decreasing this. The two stage method is not suitable for this.")
  if(first.stage.threshold == 0){
    if(stats::runif(1) < 0.5) {
      stop("first.stage.threshold is 0, consider increasing this. The two stage method is not suitable for this. \n     'None Shall Pass!' - The Black Knight (Monty Python and the Holy Grail)")
    } else
      stop("first.stage.threshold is 0, consider increasing this. The two stage method is not suitable for this. \n     'You Shall Not Pass!' - Gandalf the Grey")
  }


  if(abs(report.lowest.amount - round(report.lowest.amount)) > .Machine$double.eps^0.5 || report.lowest.amount < 1){
    warning("report.lowest.amount must be positive integer. Rounding up to positive integer")
    report.lowest.amount <- max(ceiling(report.lowest.amount), 1, na.rm = TRUE)
  }

  if(max.batchsize < 2) stop("max.batchsize must be 2 or greater.")

  if(abs(number.of.covariates - round(number.of.covariates)) > .Machine$double.eps^0.5 || number.of.covariates < 1){
    warning("number.of.covariates must be non-negative integer. Rounding up to non-negative integer")
    number.of.covariates <- max(ceiling(number.of.covariates), 0, na.rm = TRUE)
  }

  if(progress == 0) progress = FALSE

  this.call <- match.call()

  snps.are.named = !is.null(dimnames(test.covariates)[[2]])
  number.of.subjects = dim(test.covariates)[1]
  max.batchsize = dim(test.covariates)[2]

  rm(test.covariates)

  start.time <- proc.time()[3]

  if(number.of.covariates == 0){
    for(batch.num in 1:length(covariate.filepaths)){
      number.of.covariates <- number.of.covariates + dim(read.function(covariate.filepaths[batch.num]))[2]
    }
  }

  ts.output <-
    batched.secondstagecoxph(
      survival.dataset = survival.dataset,
      covariate.filepaths = covariate.filepaths,
      first.stage.threshold = first.stage.threshold,
      progress = progress,
      max.coef = max.coef,
      max.batchsize = max.batchsize,
      updatefile = updatefile,
      upper.bound.correlation = upper.bound.correlation,
      read.function = read.function,
      number.of.covariates = number.of.covariates,
      snps.are.named = snps.are.named,
      number.of.subjects = number.of.subjects
    )


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
      names.lowest.five <- matrix(c(dimnames(ts.output$second.stage.sparse.matrix)[[2]][as.vector(indices.lowest.five)]),
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
  if(progress != 0) cat("\rAnalysis completed. Runtime:", total.runtime)

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



#' Multicore method of performing the second stage
#'
#' @param survival.dataset the outcome data
#' @param first.stage.threshold the FST
#' @param progress set to 0 for no updates
#' @param max.coef maximum value of fitted weights before declared non-converged
#' @param updatefile path to that file
#' @param max.batchsize max number of covariates in one batch
#' @param upper.bound.correlation upper bound on the correlation before not checked
#' @param read.function the function used to read covariate batches
#' @param covariate.filepaths paths to the files containing the covariates
#' @param number.of.covariates combined total of covariates in all files
#' @param snps.are.named boolean, whether or not the snps will be named when read from file.
#' @param number.of.subjects Total number of subjects in all files
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
#'   passed.indices are not used here (only returned), since amount.rejections is better since
#'   possible NAs have been removed.
#'
#'   Reading from given files in combination with a maximum batchsize becomes... interesting.
#'   Firstly: note that opening a connection to a file has overhead, so we wish to
#'   avoid doing this as much as possible. This leads us to the following consideration.
#'   Since the files, and consequently the covariates within, are fixed and determined before
#'   the results of the first stage are available, we can have... suboptimal... cases. What if
#'   all covariates from one file are all rejected? Or only one for each file? What do we do then?
#'   Ideally, we would put all the rejected covariates in new files, matching the sizes of the
#'   original files (to allow for memory-usage control), but that would be non-sensical. That is,
#'   because then we would read all files individually into memory, and then write to storage again,
#'   and finally read into memory again, and only then, after all this reading and writing would
#'   the analysis begin.
#'   The solution I chose is to group subsequent files together, in such a way that the relevant,
#'   i.e. rejected in first stage, covariates from one group do not exceed the amount of covariates
#'   typically found in one file. That way the memory-usage can still be controlled by the user,
#'   by adjusting their files (which they would have done already for the first stage), while still
#'   having a reasonable performance. It is true that this leads to inefficiencies. For example, say
#'   that the first 3 batches have 40% of their covariates rejected in the first stage. Then in the
#'   second stage, we will only read the first 2 batches, resulting in 20% of the memory-capacity not
#'   being used optimally, but this is a worthwile sacrifice in my opinion.
#'   Another option would be to write or find a optimization algorithm that would rearrange the batches
#'   so the abovementioned method would waste as little memory as possible. But that algorithm would
#'   also introduce overhead, and e.g. a naive approach to this problem (try every possible permutation)
#'   would be of complexity O(n!) with n the number of batches. Maybe an efficient algorithm exists,
#'   but I'm too lazy to implement this. The number of batches should ideally be too large, so it is
#'   not that big of a deal. Anyway, that is therefore why the indices for each batch may look odd.
#'
batched.secondstagecoxph <- function(survival.dataset, covariate.filepaths, first.stage.threshold,
                                  progress = 50, max.coef = 5, updatefile = "", max.batchsize = 1000,
                                  upper.bound.correlation = 0.95,
                                  read.function = function(x) as.matrix(read.table(x)),
                                  number.of.covariates,
                                  snps.are.named = FALSE, number.of.subjects){
  first.stage.list <- batched.firststagecoxph(
    survival.dataset = survival.dataset,
    covariate.filepaths = covariate.filepaths,
    progress = progress,
    max.coef = max.coef,
    max.batchsize = max.batchsize,
    updatefile = updatefile,
    read.function = read.function,
    number.of.covariates = number.of.covariates,
    snps.are.named = snps.are.named
  )

  #unlist the result:
  first.stage.result = first.stage.list$p.values
  names(first.stage.result) = first.stage.list$names

  #revert (possible) reordering introduced in first stage
  first.stage.permutation <- order(first.stage.list$og.index)
  first.stage.result <- first.stage.result[first.stage.permutation]
  rm(first.stage.list)

  first.stage.rejections <- (first.stage.result < first.stage.threshold)
  passed.indices <- which(first.stage.rejections == TRUE)
  amount.rejections <- max(sum(first.stage.rejections, na.rm = TRUE), 0, na.rm = TRUE)

  if(amount.rejections == 0) warning("No rejections in the first stage")
  if(amount.rejections == 1) warning("Only one rejection in the first stage")
  if(amount.rejections < 2) return(list(second.stage.sparse.matrix = NA, passed.indices = passed.indices,
                                        first.stage.p.values = first.stage.result))

  no.workers <- foreach::getDoParWorkers()
  rejections.per.file.list <- list()
  for(file.count in 1:(length(covariate.filepaths)-1)){
    rejections.per.file.list = append(rejections.per.file.list, list(which(first.stage.rejections[1:max.batchsize + (file.count-1)*max.batchsize])))
  }
  rejections.per.file.list = append(rejections.per.file.list, list(which(first.stage.rejections[((length(covariate.filepaths) -
                                                                                                    1) * max.batchsize + 1):length(first.stage.rejections)])))
  rejections.per.file <- apply(
    matrix(c(first.stage.rejections, rep(FALSE, length(covariate.filepaths)*max.batchsize - length(first.stage.rejections))),
           nrow = max.batchsize, ncol = length(covariate.filepaths)),
    2, function(x) sum(x, na.rm = TRUE))

  process.file.indices <- list()
  finished = FALSE
  file.count = 1
  count.covs.in.process = 0
  files.in.this.process = numeric(0)
  for(file.count in 1:length(rejections.per.file)){
    count.covs.in.process = count.covs.in.process + rejections.per.file[file.count]
    if(count.covs.in.process <= max.batchsize){
      files.in.this.process = c(files.in.this.process, file.count)
    } else {
      process.file.indices = append(process.file.indices, list(files.in.this.process))
      files.in.this.process = file.count
      count.covs.in.process = rejections.per.file[file.count]
    }
  }
  process.file.indices = append(process.file.indices, list(files.in.this.process))
  #this is now a list, with each element of the *list* containing a row of indices. The sum of
  #amount of first stage rejections of the covariates from the files with the indices of each element,
  #do not exceed the maximum batchsize.

  #make list containing location and names of covariates
  list.location.covariates <- list(file = numeric(0), index = numeric(0), names = character(0))
  for(file.count in 1:length(covariate.filepaths)){
    covariates <- read.function(covariate.filepaths[file.count])[,rejections.per.file.list[[file.count]]]
    list.location.covariates$file = c(list.location.covariates$file, rep(file.count, rejections.per.file[[file.count]]))
    list.location.covariates$index = c(list.location.covariates$index, rejections.per.file.list[[file.count]])
    if(snps.are.named){
      list.location.covariates$names = c(list.location.covariates$names, dimnames(covariates)[[2]])
    } else {
      list.location.covariates$names = c(list.location.covariates$names, names(rejections.per.file.list[[file.count]]))
    }
  }
  #(file, index) of this list form a pseudo x-y coordinate system of the covariates for relocating them.


  start.time.second.stage <- proc.time()[3]

  #contrary to the first stage, we will not work with a matrix of indices, instead with the list of indices per file and the list of files for each process.
  first.process.index <- NULL #to suppress a note from devtools::check()
  output.matrix <- foreach::foreach(first.process.index = 1:length(process.file.indices), #no '-1' here, since the last batch must be checked with itself, contrary to covariates where this is not the case
                                    .packages = c("survival", "utils"),
                                    .export = c("prefitting.check.two", "convergence.check",
                                                "clear.current.line"),
                                    .combine = rbind) %dopar% {
    #I'm NOT gonna reindent these lines to match up with the foreach-loop :P
    first.files.process <- process.file.indices[[first.process.index]]

    first.process.covariates <- matrix(NA, nrow = number.of.subjects, ncol = 0)
    for(file.count in first.files.process){
      first.process.covariates <- cbind(first.process.covariates, read.function(covariate.filepaths[file.count])[,rejections.per.file.list[[file.count]]])
    }
    #these are now all relevant covariates. Also we have a list to remember where we can find them again.


    #first.batch.real.indices <- batch.indices.matrix[,first.process.index]
    #first.batch.real.indices <- first.batch.real.indices[!is.na(first.batch.real.indices)]

    #prepare the matrix we will return (for this first.process)
    return.matrix <- matrix(NA, nrow = dim(first.process.covariates)[2], ncol = amount.rejections)

    #first check the batch with itself, if there is more than 1 covariate in the batch
    if(dim(first.process.covariates)[2] > 1){
      for(first.local.index in 1:(dim(first.process.covariates)[2]-1)){
        for(second.local.index in (first.local.index + 1):dim(first.process.covariates)[2]){
          covariate.one <- first.process.covariates[, first.local.index]
          covariate.two <- first.process.covariates[, second.local.index]
          if(prefitting.check.two(covariate.one, covariate.two, upper.bound.correlation)){
            withCallingHandlers(
              fitted.model <-
                coxph(survival.dataset ~ covariate.one * covariate.two),
              warning = function(w) {
                #print(w)
                #print(str(w))
                if (grepl("coefficient may be infinite", w$message)) {
                  invokeRestart("muffleWarning")
                  #An error was given, which is taken into account in convergence.check
                } else if (grepl("out of iterations", w$message)) {
                  invokeRestart("muffleWarning")
                  #An error was given, which is taken into account in convergence.check
                } else {
                  message(w$message)
                }
              }
            )
            if(convergence.check(fitted.model, max.coef)){
              return.matrix[first.local.index,
                            cumsum(c(0, rejections.per.file))[first.files.process[1]] + second.local.index] <-
                summary(fitted.model)$coefficients[3,5]
              #note the discrepancy between the two subset-indices. This is intentional!
              #we need to offset the second index with the number of covariates already considered in combination
              #with the first set of covariates. This matches up with the cum.sum of the rejections per file,
              #of which we need the index that matches up to the amount of files already considered.
            }
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
    if(first.process.index != length(process.file.indices)){
      for(second.process.index in (first.process.index + 1):length(process.file.indices)){
        #I'm NOT gonna reindent these lines to match up with the foreach-loop :P
        second.files.process <- process.file.indices[[second.process.index]]

        second.process.covariates <- matrix(NA, nrow = number.of.subjects, ncol = 0)
        for(file.count in second.files.process){
          second.process.covariates <- cbind(second.process.covariates, read.function(covariate.filepaths[file.count])[,rejections.per.file.list[[file.count]]])
        }

        if(dim(second.process.covariates)[2] > 0){ #only if the second batch is not empty
          for(first.local.index in 1:dim(first.process.covariates)[2]){
            for(second.local.index in 1:dim(second.process.covariates)[2]){
              covariate.one <- first.process.covariates[, first.local.index]
              covariate.two <- second.process.covariates[, second.local.index]
              if(prefitting.check.two(covariate.one, covariate.two, upper.bound.correlation)){
                withCallingHandlers(
                  fitted.model <-
                    coxph(survival.dataset ~ covariate.one * covariate.two),
                  warning = function(w) {
                    #print(w)
                    #print(str(w))
                    if (grepl("coefficient may be infinite", w$message)) {
                      invokeRestart("muffleWarning")
                      #An error was given, which is taken into account in convergence.check
                    } else if (grepl("out of iterations", w$message)) {
                      invokeRestart("muffleWarning")
                      #An error was given, which is taken into account in convergence.check
                    } else {
                      message(w$message)
                    }
                  }
                )

                if(convergence.check(fitted.model, max.coef)){
                  return.matrix[first.local.index,
                                cumsum(c(0, rejections.per.file))[second.files.process[1]] + second.local.index] <-
                    #note the discrepancy between the two subset-indices. This is intentional!
                    #we need to offset the second index with the number of covariates already considered in combination
                    #with the first set of covariates. This matches up with the cum.sum of the rejections per file,
                    #of which we need the index that matches up to the amount of files already considered.
                    summary(fitted.model)$coefficients[3,5]
                }
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


  #output.matrix is a dense matrix, which we need to expand to a sparse matrix with proper indexing and naming
  non.na.indices <- which(!is.na(output.matrix), arr.ind = TRUE)
  #these indices are local to the output.matrix, convert these to (index.of.file-1)*max.batchsize + index.in.file
  global.index <- list.location.covariates$index + (list.location.covariates$file-1)*max.batchsize
  global.non.na.indices <- cbind(global.index[non.na.indices[,1]], global.index[non.na.indices[,2]])

  sparse.dimnames <- list(rep("NA", number.of.covariates), rep("NA", number.of.covariates))
  sparse.dimnames[[1]][global.index] = sparse.dimnames[[2]][global.index] = list.location.covariates$names

  second.stage.sparse.matrix <- Matrix::sparseMatrix(i = global.non.na.indices[,1], j = global.non.na.indices[,2],
                                                     x = as.vector(output.matrix[non.na.indices]),
                                                     triangular = TRUE, dims = rep(number.of.covariates, 2))
  dimnames(second.stage.sparse.matrix) <- sparse.dimnames


  return(list(second.stage.sparse.matrix = second.stage.sparse.matrix, passed.indices = passed.indices,
              first.stage.p.values = first.stage.result))
}



#' Perform the first stage of a GWAS in parallel
#'
#' @param survival.dataset the outcomes
#' @param progress how often progress must be reported
#' @param max.coef what the maximum coefficients may be when fitting
#' @param max.batchsize what the maximum batchsize must be.
#' @param updatefile path to file to replace terminal
#' @param covariate.filepaths paths to the files containing covariates
#' @param read.function function used to read covariates
#' @param number.of.covariates combined total of covariates in all files
#' @param snps.are.named boolean, whether or not the covariates will be named when read from file
#'
#' @return the p-values of the first stage (in random order) in a list. This list contains the
#'          p-values, the names and the original index (for reordering later).
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
batched.firststagecoxph <- function(survival.dataset, covariate.filepaths, progress = 50,
                                      max.coef = 5, max.batchsize = 1000, updatefile = "",
                                    read.function = function(x) as.matrix(read.table(x)),
                                    number.of.covariates,
                                    snps.are.named = FALSE){
  start.time.first.stage <- proc.time()[3]

  no.covariates <- number.of.covariates
  no.workers <- foreach::getDoParWorkers()

  no.processes <- ceiling(no.covariates/max.batchsize)
  batchsizes.vector <- c(rep(max.batchsize, no.processes - 1),
                         no.covariates - max.batchsize*(no.processes-1))
  matrix.indices.processes <- matrix(c(seq_len(sum(batchsizes.vector)),
                                       rep(NA, length(batchsizes.vector)*max.batchsize - sum(batchsizes.vector))),
                                     ncol = length(batchsizes.vector),
                                     nrow = max.batchsize)

  process.index <- NULL #suppressing a note from devtools::check()
  p.value.list <- foreach::foreach(process.index = seq_len(no.processes),
                                     .packages = c("survival", "utils"),
                                     .export = c("prefitting.check.one", "convergence.check",
                                                 "clear.current.line"),
                                     .combine = function(x,y) mapply(c, x, y, SIMPLIFY = FALSE)) %dopar% {
    #I will not reindent these.
    this.process.indices = (matrix.indices.processes[,process.index])[!is.na(matrix.indices.processes[,process.index])]
    if(length(this.process.indices) == 0) return(rep(NA, max.batchsize)) #if the batch is empty, return NA's

    this.batch.covariates <- read.function(covariate.filepaths[process.index])

    this.process.p.values <- rep(NA, length(this.process.indices))
    for(covariate.index in 1:length(this.process.indices)){
     this.covariate <- this.batch.covariates[,covariate.index]
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
       if (convergence.check(fitted.model, max.coef)) {
         this.process.p.values[covariate.index] <-
           summary(fitted.model)$coefficients[1, 5]
       }
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

    if(snps.are.named) {
      names(this.process.p.values) <- dimnames(this.batch.covariates)[[2]]
    } else names(this.process.p.values) <- this.process.indices
    attr(this.process.p.values, "og.index") <- this.process.indices

    #pack attributes into list:
    return.list <- list(p.values = this.process.p.values,
                        names = names(this.process.p.values),
                        og.index = this.process.indices)

    return(return.list)
                                     }


  if(progress != 0){
    clear.current.line()
    cat("\rFirst stage complete. Commencing second stage. ")
  }
  return(p.value.list)
}



