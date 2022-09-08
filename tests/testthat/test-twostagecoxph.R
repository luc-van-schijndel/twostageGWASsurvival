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

})

test_that("Incorrect first.stage.threshold input gives errors", {
  survival.dataset <- survival::Surv(1:10, rep(1, 10))
  covariate.matrix <- matrix(rep(seq(1,5), 6), 10, 3)

  expect_error(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = "string"), "numeric")
  expect_error(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 1:2), "length")
  expect_error(twostagecoxph(survival.dataset, covariate.matrix, first.stage.threshold = 2), "interval")
})

test_that("First stage function gives properly structured output", {
  survival.dataset <- survival::Surv(1:10, rep(1, 10))
  covariate.matrix <- matrix(rep(seq(1,5), 6), 10, 3)

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

  output <- twostagecoxph(survival.dataset, covariate.matrix)

  #(c(summary(survival::coxph(survival.dataset ~ covariate.matrix[,1]*covariate.matrix[,2]))$coef[3,5],
  #   summary(survival::coxph(survival.dataset ~ covariate.matrix[,1]*covariate.matrix[,3]))$coef[3,5],
  #   summary(survival::coxph(survival.dataset ~ covariate.matrix[,2]*covariate.matrix[,3]))$coef[3,5]))
  expect_equal(output[[1]], matrix(c(1, 0.29283645, 0.76400342,
                                     1, 1,          0.09613708,
                                     1, 1,          1),
                                   nrow = 3, ncol = 3, byrow = TRUE, dimnames = list(1:3, 1:3)),
               tolerance = 1e-7)
  #3x the previous (maximum 1)
  expect_equal(output[[2]], matrix(c(1, 0.8785094, 1,
                                     1, 1,         0.2884112,
                                     1, 1,         1),
                                   nrow = 3, ncol = 3, byrow = TRUE, dimnames = list(1:3, 1:3)),
               tolerance = 1e-7)
})
