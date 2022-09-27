test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

#Input error feedback==============

test_that("Wrong input objects give errors", {
  expect_error(twostagecoxph(1:2, matrix(1:4,2,2)), "Surv")
  expect_error(twostagecoxph(survival::Surv(1:2, 0:1), 1:4),  "matrix")
  expect_error(twostagecoxph(survival::Surv(1,1), matrix(1:4, 2, 2)), "No dimension")

  survival.dataset <- survival::Surv(c(5,5,3,3,2,2,2,1,1,1),
                                     c(0,0,1,1,1,1,1,1,1,1))
  covariate.matrix <- matrix(c(2,2,1,
                               2,2,1,
                               1,2,1,
                               2,1,1,
                               1,1,1,
                               1,1,1,
                               1,0,0,
                               0,1,0,
                               1,0,0,
                               0,0,0),
                             nrow = 10, ncol = 3, byrow = TRUE)
  expect_warning(twostagecoxph(survival.dataset, t(covariate.matrix), control = twostagecoxph.control(progress = 0)), "Transposed")

})

test_that("Incorrect first.stage.threshold input gives errors", {
  survival.dataset <- survival::Surv(c(5,5,3,3,2,2,2,1,1,1),
                                     c(0,0,1,1,1,1,1,1,1,1))
  covariate.matrix <- matrix(c(2,2,1,
                               2,2,1,
                               1,2,1,
                               2,1,1,
                               1,1,1,
                               1,1,1,
                               1,0,0,
                               0,1,0,
                               1,0,0,
                               0,0,0),
                             nrow = 10, ncol = 3, byrow = TRUE)

  expect_error(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = "string"), "numeric")
  expect_error(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 1:2), "length")
  expect_error(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 2), "interval")
  expect_error(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 1), "decreasing")
  expect_error(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 0), "increasing")

  expect_warning(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 1e-7, control = twostagecoxph.control(progress = 0)), "No rejections")
  expect_warning(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 0.01, control = twostagecoxph.control(progress = 0)), "one rejection") #the first stage p-values are 0.008, 0.013 and 0.015
})

test_that("Wrong input gives errors for control function", {
  expect_error(twostagecoxph.control(report.lowest.amount = c(1,2)), "integer scalar")
  expect_error(twostagecoxph.control(return.raw = c(FALSE,FALSE)), "logical scalar")
  expect_error(twostagecoxph.control(progress = c(1,2)), "integer scalar")
  expect_error(twostagecoxph.control(max.coef = c(1,2)), "integer scalar")
  expect_error(twostagecoxph.control(max.batchsize = c(1,2)), "integer scalar")
  expect_error(twostagecoxph.control(upper.bound.correlation = c(1,2)), "double scalar")
  expect_error(twostagecoxph.control(max.batchsize = 1), "2 or greater")
})

#Proper (numeric) outputs=======================

test_that("First stage function gives properly structured output", {
  survival.dataset <- survival::Surv(c(5,5,3,3,2,2,2,1,1,1),
                                     c(0,0,1,1,1,1,1,1,1,1))
  covariate.matrix <- matrix(c(2,2,1,
                               2,2,1,
                               1,2,1,
                               2,1,1,
                               1,1,1,
                               1,1,1,
                               1,0,0,
                               0,1,0,
                               1,0,0,
                               0,0,0),
                             nrow = 10, ncol = 3, byrow = TRUE)

  first.stage.output <- firststagecoxph(survival.dataset, covariate.matrix, progress = 0)

  expect_length(first.stage.output, 3)
  expect_lt(first.stage.output[1], 1)
  expect_gt(first.stage.output[1], 0)


})

test_that("Outputs are the p-values we expect", {
  survival.dataset <- survival::Surv(c(5,5,3,3,2,2,2,1,1,1),
                                     c(0,0,1,1,1,1,1,1,1,1))
  covariate.matrix <- matrix(c(2,2,1,
                               2,2,1,
                               1,2,1,
                               2,1,1,
                               1,1,1,
                               1,1,1,
                               1,0,0,
                               0,1,0,
                               1,0,0,
                               0,0,0),
                             nrow = 10, ncol = 3, byrow = TRUE)

  #(c(summary(survival::coxph(survival.dataset ~ covariate.matrix[,1]))$coef[1,5],
  #   summary(survival::coxph(survival.dataset ~ covariate.matrix[,2]))$coef[1,5],
  #   summary(survival::coxph(survival.dataset ~ covariate.matrix[,3]))$coef[1,5]))
  expected.first.stage = c(0.007898618, 0.012893229, 0.015147279)
  names(expected.first.stage) <- 1:3
  expect_equal(firststagecoxph(survival.dataset, covariate.matrix, progress = 0),
               expected.first.stage,
               tolerance = 1e-7)

  #code used to obtain constants in the following section:
  #(c(summary(survival::coxph(survival.dataset ~ covariate.matrix[,1]*covariate.matrix[,2]))$coef[3,5],
  #   summary(survival::coxph(survival.dataset ~ covariate.matrix[,1]*covariate.matrix[,3]))$coef[3,5],
  #   summary(survival::coxph(survival.dataset ~ covariate.matrix[,2]*covariate.matrix[,3]))$coef[3,5]))
  expected.second.stage <- list(second.stage.sparse.matrix = Matrix::sparseMatrix(i = c(1,1,2), j = c(2,3,3),
                                                         x = c(0.29283645, 0.76400342, 0.09613708),
                                                         triangular = TRUE),
       passed.indices = c(1,2,3),
       first.stage.p.values = c(0.007898618, 0.012893229, 0.015147279))
  names(expected.second.stage$passed.indices) = 1:3
  names(expected.second.stage$first.stage.p.values) = 1:3
  expect_equal(singlecore.twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 0.05, progress = 0),
               expected.second.stage,
               tolerance = 1e-7)

  expect_equal(twostagecoxph(survival.dataset, covariate.matrix,
                             control = twostagecoxph.control(progress = 0, return.raw = TRUE,
                                                             upper.bound.correlation = 0.95))$p.value.matrix,
               Matrix::sparseMatrix(i = c(1,1,2), j = c(2,3,3),
                                    x = c(0.29283645, 0.76400342, 0.09613708),
                                    triangular = TRUE),
               tolerance = 1e-7)

  #3x the previous (maximum 1)
  expect_equal(twostagecoxph(survival.dataset, covariate.matrix,
                             control = twostagecoxph.control(progress = 0, upper.bound.correlation = 0.95))$p.value.matrix,
               Matrix::sparseMatrix(i = c(1,1,2), j = c(2,3,3),
                                    x = c(0.8785094, 1, 0.2884112),
                                    triangular = TRUE),
               tolerance = 1e-7)

  #highest = raw, 2nd = 2x and lowest = 3x of the raw p-values
  expect_equal(twostagecoxph(survival.dataset, covariate.matrix,
                             multiple.hypotheses.correction = "hochberg",
                             control = twostagecoxph.control(progress = 0, upper.bound.correlation = 0.95))$p.value.matrix,
               Matrix::sparseMatrix(i = c(1,1,2), j = c(2,3,3),
                                    x = c(0.58567291, 0.76400342, 0.2884112),
                                    triangular = TRUE),
               tolerance = 1e-7)
})

#Multicore tests====================

test_that("(Multicore) first stage gives proper output", {
  doParallel::registerDoParallel(2)

  survival.dataset <- survival::Surv(c(5,5,3,3,2,2,2,1,1,1),
                                     c(0,0,1,1,1,1,1,1,1,1))
  covariate.matrix <- matrix(c(2,2,1,
                               2,2,1,
                               1,2,1,
                               2,1,1,
                               1,1,1,
                               1,1,1,
                               1,0,0,
                               0,1,0,
                               1,0,0,
                               0,0,0),
                             nrow = 10, ncol = 3, byrow = TRUE)

  first.stage.output <- firststagecoxph.multicore(survival.dataset, covariate.matrix, progress = 0)
  expected.output <- c(0.007898618, 0.012893229, 0.015147279)
  names(expected.output) <- c("1", "2", "3")

  expect_length(first.stage.output, 3)
  expect_lt(max(first.stage.output), 1)
  expect_gt(min(first.stage.output), 0)
  expect_equal(first.stage.output[c("1", "2", "3")], expected.output,
               tolerance = 1e-7)

  #This tests the case where one batch is completely empty:
  first.stage.output <- firststagecoxph.multicore(survival.dataset, covariate.matrix,
                                                  progress = 0, max.batchsize = 1)
  expected.output <- c(0.007898618, 0.012893229, 0.015147279, NA)
  names(expected.output) <- c("1", "2", "3", "")

  expect_length(first.stage.output, 4)
  expect_lt(max(first.stage.output, na.rm = TRUE), 1)
  expect_gt(min(first.stage.output, na.rm = TRUE), 0)
  expect_equal(first.stage.output[c("1", "2", "3")], expected.output[c("1", "2", "3")],
               tolerance = 1e-7) #we excluded subsetting with "", since the name attribute then gets lost.

  #second stage multicore function---------------

  second.stage.output <- multicore.twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 0.05, progress = 0)
  second.stage.output.singlecore <- singlecore.twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 0.05, progress = 0)
  expect_equal(second.stage.output, second.stage.output.singlecore)
})

#optimal batch configuration tests================

test_that("Optimal batch configuration calculations are right", {
  #For first stage (normal):

  #When restrictions are not a problem:
  expect_equal(optimal.batch.configuration(no.covariates = 10, max.batchsize = 1000, no.workers = 2),
               list(optimal.batchsize = 5, optimal.no.batches = 2, last.batchsizes = c(5,5)))
  expect_equal(optimal.batch.configuration(no.covariates = 20, max.batchsize = 1000, no.workers = 3),
               list(optimal.batchsize = 7, optimal.no.batches = 3, last.batchsizes = c(7,7,6)))
  expect_equal(optimal.batch.configuration(no.covariates = 30, max.batchsize = 1000, no.workers = 32),
               list(optimal.batchsize = 1, optimal.no.batches = 32, last.batchsizes = c(rep(1,30), rep(0, 2))))

  #When restrictions are a problem:
  expect_equal(optimal.batch.configuration(no.covariates = 50, max.batchsize = 10, no.workers = 2),
               list(optimal.batchsize = 9, optimal.no.batches = 6, last.batchsizes = c(7,7))) #c(9, 9, 9, 9, 7, 7) is better than c(10, 10, 10, 10, 10, 0)
  expect_equal(optimal.batch.configuration(no.covariates = 100, max.batchsize = 20, no.workers = 3),
               list(optimal.batchsize = 17, optimal.no.batches = 6, last.batchsizes = c(17, 17, 15)))
  expect_equal(optimal.batch.configuration(no.covariates = 300, max.batchsize = 10, no.workers = 16),
               list(optimal.batchsize = 10, optimal.no.batches = 32, last.batchsizes = c(rep(9, 15), 5)))

  #When edge-cases:
  expect_equal(optimal.batch.configuration(no.covariates = 50, max.batchsize = 2, no.workers = 4),
               list(optimal.batchsize = 2, optimal.no.batches = 28, last.batchsizes = c(1, 1, 0, 0)))
  expect_equal(optimal.batch.configuration(no.covariates = 100, max.batchsize = 2, no.workers = 3),
               list(optimal.batchsize = 2, optimal.no.batches = 51, last.batchsizes = c(2, 2, 0)))
  expect_equal(optimal.batch.configuration(no.covariates = 300, max.batchsize = 10, no.workers = 32),
               list(optimal.batchsize = 10, optimal.no.batches = 32, last.batchsizes = c(rep(10, 30), rep(0, 2))))


  #Using it in the second stage ---------
  #we divide max.batchsize and no.covariates by 2 (rounding down). max.batchsize since we need to compare (mostly) 2 batches,
  #so we halve this limit to accomodate, and the number of covariates, since we split the process
  #up in. This allows for similar finish-times using first-in, last-out (see documentation of multicore.twostagecoxph)

  #When restrictions are not a problem:
  expect_equal(optimal.batch.configuration(no.covariates = 10/2, max.batchsize = 1000/2, no.workers = 2),
               list(optimal.batchsize = 3, optimal.no.batches = 2, last.batchsizes = c(3,2))) #1:3x(1:3, 4:5, 6:8, 9:10) + 4:5x(4:5, 6:8, 9:10) + 6:8x(6:8, 9:10) + 9:10x9:10
  expect_equal(optimal.batch.configuration(no.covariates = 20/2, max.batchsize = 1000/2, no.workers = 3),
               list(optimal.batchsize = 4, optimal.no.batches = 3, last.batchsizes = c(4,4,2)))
  expect_equal(optimal.batch.configuration(no.covariates = 30/2, max.batchsize = 1000/2, no.workers = 32),
               list(optimal.batchsize = 1, optimal.no.batches = 32, last.batchsizes = c(rep(1,15), rep(0, 17))))

  #When restrictions are a problem:
  expect_equal(optimal.batch.configuration(no.covariates = 50/2, max.batchsize = 10/2, no.workers = 2),
               list(optimal.batchsize = 5, optimal.no.batches = 6, last.batchsizes = c(3,2))) #c(9, 9, 9, 9, 7, 7) is better than c(10, 10, 10, 10, 10, 0)
  expect_equal(optimal.batch.configuration(no.covariates = 100/2, max.batchsize = 20/2, no.workers = 3),
               list(optimal.batchsize = 9, optimal.no.batches = 6, last.batchsizes = c(8,8,7)))#same as above, batches of 9 distribute the load better.
  expect_equal(optimal.batch.configuration(no.covariates = 300/2, max.batchsize = 10/2, no.workers = 16),
               list(optimal.batchsize = 5, optimal.no.batches = 32, last.batchsizes = c(rep(5, 14), rep(0, 2))))

  #When edge-cases:
  expect_equal(optimal.batch.configuration(no.covariates = 50/2, max.batchsize = 2/2, no.workers = 4),
               list(optimal.batchsize = 1, optimal.no.batches = 28, last.batchsizes = c(1, 0, 0, 0)))
  expect_equal(optimal.batch.configuration(no.covariates = 100/2, max.batchsize = 2/2, no.workers = 3),
               list(optimal.batchsize = 1, optimal.no.batches = 51, last.batchsizes = c(1, 1, 0)))
  expect_equal(optimal.batch.configuration(no.covariates = 300/2, max.batchsize = 10/2, no.workers = 32),
               list(optimal.batchsize = 5, optimal.no.batches = 32, last.batchsizes = c(rep(5, 30), rep(0, 2))))
})
