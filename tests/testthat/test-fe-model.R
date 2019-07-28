context("FE Model: Simple Fixed Effects")

source("plain-data.R")

test_that("Demean should be correct w.r.t. felm.", {
  library(lfe)

  ## inds     <- create.indicators(raw.inds)
  actual   <- unname(fastplm(ind = raw.inds, y = y, x = x)$coefficients)
  expected <- unname(felm(y ~ x | raw.inds[,1] + raw.inds[,2] + raw.inds[,3])$coefficients)

  expect_equal(actual, expected)
})

test_that("Demean should be independent of concrete values of indicators.", {
  ## inds1 <- create.indicators(raw.inds)
  ## inds2 <- create.indicators(new.inds)

  expected <- unname(fastplm(ind = raw.inds, y = y, x = x)$coefficients)
  actual   <- unname(fastplm(ind = new.inds, y = y, x = x)$coefficients)

  expect_equal(actual, expected)
})

test_that("Estimated fixed effects should have correct row names.", {
  ## inds  <- create.indicators(new.inds)
  model <- fastplm(ind = new.inds, y = y, x = x)
  inds <- model$inds

  sfe.coefs <- lapply(model$sfe.coefs, function(coefs) as.list(coefs[, 1]))
  effs      <- Reduce("+", lapply(1 : length(sfe.coefs), function(col) {
    sapply(1 : nrow(y), function(row) {
      (sfe.coefs[[col]])[[inds$inds[row, col]]]
    })
  }))

  fitted.values <- x %*% model$coefficients + model$intercept + effs

  expect_equal(fitted.values, model$fitted.values)
})

test_that("Demean without X should yield fixed effects.", {
  ## inds  <- create.indicators(new.inds)
  model <- fastplm(y = effs, ind = new.inds)

  for (i in 1 : 3)
    expect_equal(unname(model$sfe.coefs[[i]]),
                 unname(as.matrix(sfe.coefs[,i] - mean(sfe.coefs[,i]))))

  expect_equal(model$intercept, sum(sapply(1 : 3, function(i) mean(sfe.coefs[,i]))))
})

test_that("Parallel demean should give identical results.", {
  ## inds   <- create.indicators(new.inds)
  models <- lapply(1 : 8, function(core.num) {
    model <- fastplm(ind = raw.inds, y = y, x = x, core.num = core.num)
    # "fe" contains a pointer.
    model$fe <- NULL
    # uid is a random number
    model$inds$uid <- NULL
    model
  })

  baseline <- models[[1]]
  lapply(2 : 8, function(i) expect_identical(models[[i]], baseline))
})

test_that("Prediction with original input should yield fitted values.", {
  ## inds     <- create.indicators(new.inds)
  model    <- fastplm(ind = new.inds, y = y, x = x)
  ## sub.inds <- create.subindicators(new.inds, model)
  actual   <- predict(model, x = x, ind = new.inds)

  expect_equal(actual, model$fitted.values)
})

test_that("Prediction should ignore NA row(s).", {
  rows     <- sample(N, 5)

  ## inds     <- create.indicators(new.inds)
  model    <- fastplm(ind = new.inds, y = y, x = x)
  ## sub.inds <- create.subindicators(new.inds, model)
  expected <- predict(model, x = x, ind = new.inds)

  expected[rows, ] <- NA

  act.inds <- new.inds
  act.inds[rows, ] <- NA

  ## sub.inds <- suppressWarnings(create.subindicators(act.inds, model))
  actual <- predict(model, x = x, ind = act.inds)

  expect_identical(actual, expected)
})

context("FE Model: Complex Fixed Effects")

source("complex-effects.R")

## inds <- create.indicators(raw.inds)

## cfe1 <- create.complex.effect(inds, 1, 2, t(as.matrix(inf2)))
## cfe2 <- create.complex.effect(inds, 2, 1, t(as.matrix(inf1)))

## fe   <- create.fixed.effects(inds, cfes = list(cfe1, cfe2))

test_that("Demean should be correct w.r.t. lm.", {
  model.us <- fastplm(ind = cbind(mapped.inf(1), mapped.inf(2)), y = y, x = x, 
                      cfe = list(c(1,2),c(2,1)), PCA = FALSE)
  model.lm <- lm(y ~ x + factor(raw.inds[, 1]) * mapped.inf(2)
                       + factor(raw.inds[, 2]) * mapped.inf(1))

  actual   <- unname(model.us$coefficients)
  expected <- as.matrix(unname(model.lm$coefficients[2 : 6]))

  expect_equal(actual, expected)
})

test_that("Prediction with original input should yield fitted values.", {
  model    <- fastplm(ind = cbind(mapped.inf(1),mapped.inf(2)), y = y, x = x, 
                      cfe = list(c(1,2),c(2,1)), PCA = FALSE)
  ## sub.inds <- create.subindicators(raw.inds, model)
  actual   <- predict(model, x = x, ind = cbind(mapped.inf(1),mapped.inf(2)))

  expect_equal(actual, model$fitted.values)
})
