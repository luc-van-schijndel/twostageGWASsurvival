# twostageGWASsurvival
An R package that implements a two stage statistical method for detecting epistasis in GWAS for survival data which follows the Cox proportional hazards model. 

Installation of the package can be done by using the `install_github` function from the `devtools` package:

```
devtools::install_github("luc-van-schijndel/twostageGWASsurvival", build_vignettes = TRUE)
library(twostageGWASsurvival)
```

For a guide on its use, see the vignette html file in the vignettes folder, or use the command `vignette('twostageGWASsurvival')` after installation. 
The package has functionality for parallel processing, as well as for the user to provide the covariates in separate files to alleviate memory issues. 

The theory for the validity of this statistical test is described in a paper by Jonker et al., which can be found [here](404 not found) (paper to be published as of 09-06-2024).

