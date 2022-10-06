test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

#batched_twostagecoxph() uses the same tests as twostagecoxph() to check for proper input, so
#these tests will not be added here.

test_that("batched variant gives same output as normal variant", {
  skip_on_cran()

  doParallel::registerDoParallel(2)
  survival.dataset <- example_survival_data
  covariate.matrix <- example_snp_data
  output.direct <- twostagecoxph(survival.dataset, covariate.matrix,
                                 control = twostagecoxph.control(progress = 0))

  number.of.covs <- dim(covariate.matrix)[2]
  number.of.files <- 6
  temp.snpfile.paths <- tempfile(rep("snpfile", number.of.files), tmpdir = tempdir(check = TRUE), fileext = ".txt")
  indices.matrix <- matrix(c(seq_len(number.of.covs), rep(NA, ceiling(number.of.covs/number.of.files)*number.of.files - number.of.covs)),
                           ncol = number.of.files)
  for(file.num in 1:number.of.files){
    indices <- indices.matrix[,file.num]
    write.table(covariate.matrix[,indices[!is.na(indices)]], file = temp.snpfile.paths[file.num])
  }


  output.batched <- batched.twostagecoxph(survival.dataset, temp.snpfile.paths,
                               number.of.covariates = number.of.covs, control =
                                 twostagecoxph.control(progress = 0))

  #these objects are inherently different, since the call and the runtime differ.
  #so we check only the relevant parts
  expect_equal(output.batched$fst, output.direct$fst)
  expect_equal(output.batched$marginal.significant, output.direct$marginal.significant)
  expect_equal(output.batched$first.stage, output.direct$first.stage)

  expect_equal(output.batched$p.value.matrix, output.direct$p.value.matrix)

  expect_equal(output.batched$most.significant.results, output.direct$most.significant.results)

})
