#*batch making for multicore, either by file or by subsetting (D)
#*
#*check which calls to functions are made in double-for loops, and remove :: operators by importFrom-ing (D)
#*
#*allow for the option of never having the full covariate matrix loaded in memory. Instead reading
#*bit by bit from files. (pre-work done) (D)
#*
#*Additional return info (D, I think? No more can be reasonably added?)
#*
#*large dataset. Make good example from this. (D)
#*
#*Non-convergence checks (partially done?)
#*
#*Edge-case handling (partially done?)
#*
#*Coxph properly calling
#*Requires a lot of manual setting of parameters. Not worth the effort I think.
#*
#*Document datasets (D)
#*
#*It is possible for the second stage batches to be non-ordered returned. Maybe permutate back,
#*or require the ordered-finishing of processes? (D)
#*
#*
#*Parallel computing (D)
#*
#*Consider reparameterization of double for-loop in second stage: iterate over "relevant.indices"?
#*It seems as if the double for-loop is faster. (D)
#*
#*Consider removing the corrected spars (D)
#*
#*Give print method proper lay-out and structure (and documentation) (D)
#*
#*Give summary method proper lay-out and structure (and documentation). Display (relevant) top part of matrix (passed)
#*
#*implement ETA returning (D)
#*
#*Update documentation to match standard style: "data type"  "thing it determines in method" (D)
#*
#*Check version requirements for Depends packages. Also check if they actually need to be there. (D? I think?)
#*
#*#Checking the package should take as little CPU time as possible, as the CRAN check farm is a very limited resource and there are thousands of packages. Long-running tests and vignette code can be made optional for checking, but do ensure that the checks that are left do exercise all the features of the package.
#If running a package uses multiple threads/cores it must never use more than two simultaneously: the check farm is a shared resource and will typically be running many checks simultaneously.
# vignette niet checken

# examples should run no more than a few seconds each. Batched.twostagecoxph aanpassen. (D)

# check that everything works on other operating systems. R CMD check?

# R CMD check --as-cran ON THE TARBALL to be submitted. with current version of R-devel
# https://cran.r-project.org/web/packages/policies.html

# add ... to functions to be passed to coxph() (D)

# system.time is unreliable for child processes on non-unix likes. Documentations. (D)

#try Rprof(). See ?utils::Rprof

#Url terug doen in DESCRIPTION : URL: https://github.com/luckylluck2/twostageGWASsurvival
#                                BugReports: https://github.com/luckylluck2/twostageGWASsurvival/issues (D)

#check temporary files management. Is dit goed gedaan? unlink() files na example?

#un-comment lines of vignette

#ld-plot in examples (?) en vignette. Packages in Suggests (?)

#tests voor doParallel niet geinstalleerd-checks toevoegen.

# check usage of optimal batch conf to update output

# parallel first stage seems to allocate weird covariates

# add timestamps to logs (D)
