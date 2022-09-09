test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

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
  expect_warning(twostagecoxph(survival.dataset, t(covariate.matrix)), "Transposed")
  expect_warning(twostagecoxph(survival.dataset, covariate.matrix, report.lowest.amount = 9.5), "integer")

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

  expect_error(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 1e-7), "No rejections")
  expect_error(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 0.01), "one rejection") #the first stage p-values are 0.008, 0.013 and 0.015
})

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

  expect_length(firststagecoxph(survival.dataset, covariate.matrix), 3)
  expect_lt(firststagecoxph(survival.dataset, covariate.matrix)[1], 1)
  expect_gt(firststagecoxph(survival.dataset, covariate.matrix)[1], 0)
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
  expect_equal(firststagecoxph(survival.dataset, covariate.matrix),
               c(0.007898618, 0.012893229, 0.015147279),
               tolerance = 1e-7)


  #(c(summary(survival::coxph(survival.dataset ~ covariate.matrix[,1]*covariate.matrix[,2]))$coef[3,5],
  #   summary(survival::coxph(survival.dataset ~ covariate.matrix[,1]*covariate.matrix[,3]))$coef[3,5],
  #   summary(survival::coxph(survival.dataset ~ covariate.matrix[,2]*covariate.matrix[,3]))$coef[3,5]))
  expect_equal(twostagecoxph(survival.dataset, covariate.matrix, return.raw = TRUE)$p.value.matrix,
               Matrix::sparseMatrix(i = c(1,1,2), j = c(2,3,3),
                                    x = c(0.29283645, 0.76400342, 0.09613708),
                                    triangular = TRUE),
               tolerance = 1e-7)

  #3x the previous (maximum 1)
  expect_equal(twostagecoxph(survival.dataset, covariate.matrix)$p.value.matrix,
               Matrix::sparseMatrix(i = c(1,1,2), j = c(2,3,3),
                                    x = c(0.8785094, 1, 0.2884112),
                                    triangular = TRUE),
               tolerance = 1e-7)

  #highest = raw, 2nd = 2x and lowest = 3x of the raw p-values
  expect_equal(twostagecoxph(survival.dataset, covariate.matrix, multiple.hypotheses.correction = "hochberg")$p.value.matrix,
               Matrix::sparseMatrix(i = c(1,1,2), j = c(2,3,3),
                                    x = c(0.58567291, 0.76400342, 0.2884112),
                                    triangular = TRUE),
               tolerance = 1e-7)
})
