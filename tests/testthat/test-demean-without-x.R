context("Demean without X")

source("make-small-data.R")

test_that("Demean without X should yield fixed effects.", {
  list2env(make.small.data(), environment())
  with.effects <- function(j) sapply(1 : N, function(i) effs[inds[i, j], j])

  y <- rowSums(sapply(1 : 3, with.effects))
  result <- solveFE(as.matrix(y), inds)

  sapply(1 : 3, function(i) {
    expect_equal(result$FEcoefs[[i]],
                 as.matrix(effs[,i] - mean(effs[,i])))
  })

  expect_equal(result$intercept, sum(sapply(1 : 3, function(i) mean(effs[,i]))))
})
