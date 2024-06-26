% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twostagecoxph.R
\name{twostagecoxph.control}
\alias{twostagecoxph.control}
\title{Ancillary arguments for controlling the operation of the two stage method using twostagecoxph}
\usage{
twostagecoxph.control(
  report.lowest.amount = 5,
  return.raw = FALSE,
  progress = 1000,
  max.coef = 5,
  max.batchsize = 1000,
  lower.bound.variance = 0.1,
  upper.bound.correlation = 0.9
)
}
\arguments{
\item{report.lowest.amount}{integer, default 5; denotes how many of the most significant interactions
the function should report. Altering this value does
not affect the maximum amount reported by the \code{print.twostageGWAS} function.}

\item{return.raw}{logical, default FALSE; whether or not the output should contain the raw p-values or
the multiple hypotheses corrected p-values.}

\item{progress}{numeric, default 1000;  how many iterations should pass silently until an update is given
about runtime and progress until completion of stages. Set to 0 for no output.}

\item{max.coef}{numeric, default 5; maximum value for all coefficients in the fitted models. If any
are larger than this (in absolute value), then the model rejects them and ignored in
further analysis. Can be used to exclude unrealistically large values.}

\item{max.batchsize}{maximum size of one batch, default 1000; for parallel computing,
should be lowered if issues arise concerning memory. Can be raised
to slightly increase performance}

\item{lower.bound.variance}{numeric scalar in the interval \code{[0,1)}, default 0.1; the lower bound on the variance of a covariate
in the first stage. If the variance is lower, the model using this covariate is not fitted and the covariate is
rejected by default.}

\item{upper.bound.correlation}{numeric scalar in the interval \code{(0,1]}, default 0.9; the upper bound on the correlation between two
covariates and their interaction. If exceeded (in absolute value), the model using these two correlated variables
is not fitted and the interaction between these is rejected by default.}
}
\value{
a list containing the values of each of the above constants.
}
\description{
This is used to set various numeric parameters controlling the operation of the
               two stage method performed by \code{twostagecoxph}. Typically it would only be used
               in a call to \code{twostagecoxph}.
}
