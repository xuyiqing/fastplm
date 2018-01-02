context("Predict Fixed Effects")

source("make-small-data.R")

list2env(make.small.data(), environment())

model.1.to.n     <- solveFE(cbind(y, x), inds, estimateFE = TRUE)
predicted.1.to.n <- predictFE(model.1.to.n, x, inds)

inds <- change.group.indicators(inds)

model.chars      <- solveFE(cbind(y, x), inds, estimateFE = TRUE)
predicted.chars  <- predictFE(model.chars, x, inds)

test_that("Fitted values should be predicted with originals X where group indicators start from 1.", {
  expect_equal(model.1.to.n$fitted.values, predicted.1.to.n)
})

test_that("Fitted values should be predicted with originals X where group indicators are arbitrary.", {
  expect_equal(model.chars$fitted.values, predicted.chars)
})

test_that("predictFE should give equal results for isomorphic group indicators.", {
  expect_equal(predicted.1.to.n, predicted.chars)
})
