load("data/example.snp.data.rda") #"example.snp.data"
example.snp.data <- read.table("data-raw/s1_gendat.txt", header = TRUE)[,1:3000]


(epistasis.snps <- c(49, 485))
subjects.subset <- 1:200
covariates.subset <- 1:500
check.cor <- c(epistasis.snps-1, epistasis.snps, epistasis.snps+1)
stats::cor(example.snp.data[,check.cor])
set.seed(37 * 42)
gamma0 <- 0
gamma3 <- 0.5
time.to.event <-
  stats::rweibull(dim(example.snp.data)[1], shape = 2, scale = 1) *
  exp(gamma0 + gamma3 * example.snp.data[, epistasis.snps[1]] * example.snp.data[, epistasis.snps[2]])

censoring.time <-
  runif(dim(example.snp.data)[1],
        min = 0,
        max = quantile(time.to.event, probs = 0.75))

example.survival.data <-
  survival::Surv(pmin(time.to.event, censoring.time),
                 time.to.event <= censoring.time)

cox.model <- survival::coxph(example.survival.data[subjects.subset,] ~ example.snp.data[subjects.subset, epistasis.snps[1]] * example.snp.data[subjects.subset, epistasis.snps[2]])
print(summary(cox.model)$coef[3,5])

example_survival_data <- example.survival.data[subjects.subset,]
example_snp_data <- as.matrix(example.snp.data[subjects.subset, covariates.subset])
attr(example_snp_data, "dimnames") = list(subject.numbers = attr(example.snp.data, "dimnames")[[1]],
                                          snp.names = attr(example.snp.data, "dimnames")[[2]])

foo <- twostagecoxph(example_survival_data, example_snp_data, progress = 50, first.stage.threshold = 0.00001)
print(foo)

usethis::use_data(example_survival_data, compress = "bzip2", overwrite = TRUE)
usethis::use_data(example_snp_data, compress = "bzip2", overwrite = TRUE)

