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
#*    #return(ceiling(covs/(ceiling(ceiling(covs/min(ceiling(covs/cores), max))/cores)*cores)))
#*    }
# print("")
#
# opt.batch <- function(covs, max, cores = 2){
#   print(inter.batchsize <- min(ceiling(covs/cores), max))         #calculate initial ideal batchsize
#   print(inter.no.batches<- ceiling(covs/inter.batchsize))         #calculate how many batches this will give
#   print(final.iters     <- ceiling(inter.no.batches/cores))       #the number of iterations it will take
#   print(final.no.batches<- final.iters*cores)                     #round up to multiple of cores
#   print(final.batchsize <- ceiling(covs/final.no.batches))        #then find the corresponding batchsize
#   print(final.batchsize*final.no.batches >= covs)
#   print(final.batchsize*(final.iters - 1) < covs)
#
#   remaining.covs = covs
#   if(final.iters > 1){
#     for(i in 1:(final.iters-1)){
#       cat("Processed:", rep(final.batchsize, cores), "\n")
#       remaining.covs = remaining.covs - cores*final.batchsize
#       cat("Remaining:", remaining.covs, "\n")
#     }
#   }
#
#   print(last.iter.covs <- covs - final.batchsize*(final.iters-1)*cores)
#   last.batchsize <- ceiling(last.iter.covs/cores)
#   last.no.full.batches <- ceiling(last.iter.covs/last.batchsize - 1)
#   last.partial.batch <- last.iter.covs - last.batchsize*last.no.full.batches
#
#   cat("Processed:", c(rep(last.batchsize, last.no.full.batches), last.partial.batch), "\n")
#   remaining.covs = remaining.covs - last.batchsize*last.no.full.batches - last.partial.batch
#   cat("Remaining:", remaining.covs, "\n")
#   return(list(optimal.batchsize = final.batchsize, optimal.no.batches = final.no.batches,
#               last.batchsizes = c(rep(last.batchsize, last.no.full.batches), last.partial.batch)))
# }




#extract stuff and make matrix, with corresponding dimnames

