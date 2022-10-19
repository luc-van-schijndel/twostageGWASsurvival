## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'R6'

  R6 is a build-time dependency.

All code was hand-written by myself, so no credits needs to be given?


tests call doParallel::registerDoParallel(2), and the code itself calls foreach::getDoParWorkers(). I assume that the call to getDoParWorkers does not interfere with the parallel testing of multiple packages of CRAN. 


In one example, files are made and read using paths obtained from tempfile(). They are not removed after the example ends.
