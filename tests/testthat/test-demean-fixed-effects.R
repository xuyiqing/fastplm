library(lfe)
library(fastplm)
context("Demean Fixed Effects")

SEED  <- 19260817
N     <- 2000
LEVEL <- 50

make.data <- function () {
  set.seed(SEED)

  x <- matrix(rnorm(N * 5, 3), N, 5)
  e <- matrix(rnorm(N, 1), N, 1)

  inds <- matrix(sample(LEVEL, N * 3, replace = TRUE), N, 3)
  effs <- matrix(runif(LEVEL * 3), LEVEL, 3)

  with.effects <- function(j) sapply(1 : N, function(i) effs[inds[i, j], j])

  y <- 7 * x[,1] + 3 * x[,2] + 2 * x[,3] + 5 * x[,4] + 8 * x[,5] + 5 + e +
       with.effects(1) + with.effects(2) + with.effects(3)

  mget(c("y", "x", "e", "inds", "effs"), environment())
}

list2env(make.data(), environment())

result.lfe      <- felm(y ~ x | inds[,1] + inds[,2] + inds[,3])
result.fastplm1 <- solveFE(cbind(y, x), inds)

change.group.values <- function (inds) {
  prefixes <- c("alpha", "beta", "gamma")
  inds <- sapply(1 : 3, function(i) {
    rename <- function(j) paste(prefixes[[i]], inds[j, i], sep = ".")
    sapply(1 : N, rename)
  })
}

inds <- change.group.values(inds)

result.fastplm2 <- solveFE(cbind(y, x), inds)

test_that("Demean fixed effects where group values start from 1.", {
  expect_equal(unname(result.lfe$coefficients),
               unname(result.fastplm1$coefficients))
})

test_that("Demean fixed effects where group values are arbitrary.", {
  expect_identical(unname(result.fastplm1$coefficients),
                   unname(result.fastplm2$coefficients))
})
