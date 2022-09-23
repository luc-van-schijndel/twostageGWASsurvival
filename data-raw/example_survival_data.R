## code to prepare `example_survival_data` dataset goes here

load("data/example_snp_data.rda") #"example.snp.data"
example.snp.data <- read.table("data-raw/s1_gendat.txt", header = TRUE)[,1:3000]

seed.var = 1
seed.var = seed.var + 1
set.seed(seed.var)

<<<<<<< HEAD

(epistasis.snps <- c(2074, 2306))#sample(2000:2500, 2))
=======
(epistasis.snps <- sample(2000:2500, 2))
>>>>>>> 6648e4a67c3cf5404323bb0e9d48417de3dd50aa
subjects.subset <- 1:200
covariates.subset <- 2000:2500
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
example_snp_data <- as.matrix(example.snp.data[subjects.subset, covariates.subset], dimnames =
                                list(subject.numbers = 1:dim(example.snp.data)[1],
                                     snp.names = attr(example.snp.data, "names")[[2]]))

foo <- twostagecoxph(example_survival_data, example_snp_data, progress = 50, first.stage.threshold = 0.0001)
print(foo)

usethis::use_data(example_survival_data, compress = "bzip2", overwrite = TRUE)
usethis::use_data(example_snp_data, compress = "bzip2", overwrite = TRUE)

