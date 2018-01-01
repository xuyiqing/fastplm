library(lfe)
library(fastplm)
context("Demean Fixed Effects")

source("make-small-data.R")

list2env(make.small.data(), environment())

result.lfe      <- felm(y ~ x | inds[,1] + inds[,2] + inds[,3])
result.fastplm1 <- solveFE(cbind(y, x), inds)

inds <- change.group.indicators(inds)

result.fastplm2 <- solveFE(cbind(y, x), inds)

test_that("Demean fixed effects where group indicators start from 1.", {
  expect_equal(unname(result.lfe$coefficients),
               unname(result.fastplm1$coefficients))
})

test_that("Demean fixed effects where group indicators are arbitrary.", {
  expect_identical(unname(result.fastplm1$coefficients),
                   unname(result.fastplm2$coefficients))
})
