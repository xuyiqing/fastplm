context("Predict Fixed Effects with Unknown Indicators")

source("make-small-data.R", local = TRUE)

N     <- 50
LEVEL <- 5

list2env(make.small.data(), environment())
inds      <- change.group.indicators(inds)
model     <- solveFE(cbind(y, x), inds, estimateFE = TRUE)
predicted <- predictFE(model, x, inds)

lock.out.and.compare <- function (row) {
  col        <- (row - 1) %% 3 + 1
  inds_      <- replace(inds, (col - 1) * N + row, "unknown")
  predicted_ <- expect_warning(predictFE(model, x, inds_),
                  warning.predict.fixed.effects.with.unknown.indicators(row, col, inds_),
                  fixed = TRUE)
  expected   <- replace(predicted, row, NA)
  expect_identical(predicted_, expected)
}

test_that("Observation with unknown indicator should be predicted as N/A while not affecting others.", {
  sapply(1 : N, lock.out.and.compare)
})
