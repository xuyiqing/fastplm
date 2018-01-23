context("Estimated Fixed Effects Rownames")

source("make-small-data.R")

list2env(make.small.data(), environment())

inds <- change.group.indicators(inds)

result <- solve.fixed.effects(cbind(y, x), inds)

test_that("Estimated fixed effects interpreted using row names should yield correct fitted values.", {
  effs <- lapply(result$FEcoefs, function(coefs) as.list(coefs[,1]))

  fixed.effects <- Reduce("+", lapply(1 : ncol(inds), function(col) {
    sapply(1 : nrow(inds), function(row) {
      (effs[[col]])[[inds[row, col]]]
    })
  }))

  expect_equal(x %*% result$coefficients + result$intercept + fixed.effects,
               result$fitted.values)
})
