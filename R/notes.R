#*
#* Matrix::sparceMatrix(i = "de rij indices van de elementen",
#*                      p = "hoeveel elementen er ZIJN GEWEEST na iedere kolom (incl. de 0-kolom die niet bestaat)"
#*                      )
#*
#*
#*(test.indices <- combn(c(1,2,3),2))
#*M <- Matrix::sparseMatrix(i = test.indices[1,], j = test.indices[2,], x = 0.01)
#*str(M)
#*M
#*p.adjust(M@x)
#*
#*attr(lowest, "interaction name") = c(apply(Matrix::which(A < 20 & A > 0, arr.ind = TRUE),1, paste0(., collapse = "x")))
#*for(i in 1:10){
#*cat("\r", i, proc.time())
#*utils::flush.console()
#*Sys.sleep(1)
#*}
#*dim = 10
#*print("test")
#*for(i in 1:(dim-1)){
#*  for(j in (i+1):dim){
#*      #print(c(dim*(dim-1)/2 - (dim-i)*(dim-i-1)/2 - (dim-j), dim-i, dim-j))
#*          print(c(dim*(i-1) + j - i*(i+1)/2  )/(dim*(dim-1)/2))
#*      }
#*    }
#*    cbz <- function(max, covs, cores = 2){
#*    (intermediate.bz <- min(ceiling(covs/cores), max))
#*    (intermediate.nb <- ceiling(covs/intermediate.bz))
#*    (final.nb <- ceiling(intermediate.nb/cores)*cores)
#*    return(ceiling(covs/final.nb))
#*    #return(ceiling(covs/ceiling(ceiling(covs/min(ceiling(covs/cores), max))/cores)*cores))
#*    }
