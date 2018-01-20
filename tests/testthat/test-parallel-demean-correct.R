context("Parallel Demean Fixed Effects")

source("make-small-data.R")

list2env(make.small.data(), environment())

max.core.num <- 8

results <- lapply(1 : max.core.num,
  function(core.num) solveFE(cbind(y, x), inds, core.num))

test_that("Parallel demean should give identical results independent of cores used.", {
  baseline <- results[[1]]

  lapply(results, function(result) {
    expect_identical(result, baseline)
  })
})
