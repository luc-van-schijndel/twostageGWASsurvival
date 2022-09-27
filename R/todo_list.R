#*batch making for multicore, either by file or by subsetting
#*
#*check which calls to functions are made in double-for loops, and remove :: operators by importFrom-ing
#*
#*allow for the option of never having the full covariate matrix loaded in memory. Instead reading
#*bit by bit from files.
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
#*Document datasets
#*
#*Parallel computing (D)
#*
#*Consider reparameterization of double for-loop in second stage: iterate over "relevant.indices"?
#*It seems as if the double for-loop is faster. (D)
#*
#*Consider removing the corrected spars (D)
#*
#*Give print method proper lay-out and structure (and documentation)
#*
#*Give summary method proper lay-out and structure (and documentation). Display (relevant) top part of matrix
#*
#*implement ETA returning (D)
#*
#*Update documentation to match standard style: "data type"  "thing it determines in method"
#*
#*Check version requirements for Depends packages. Also check if they actually need to be there.
