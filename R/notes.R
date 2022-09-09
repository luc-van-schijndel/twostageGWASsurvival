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
