context("Component Analysis")

SEED <- 19260817
set.seed(SEED)

company.fe <- as.list(rnorm(4, 10))
names(company.fe) <- c("facebook", "microsoft", "amazon", "google")

person.fe <- as.list(rnorm(4, 10))
names(person.fe) <- c("alice", "bob", "cathy", "david")

inds <- rbind(c("facebook", "alice"),
              c("facebook", "bob"),
              c("microsoft", "bob"),
              c("amazon", "cathy"),
              c("amazon", "david"),
              c("google", "david"))

make.small.data <- function() {
  N <- 6
  x <- matrix(rnorm(N * 2, 3), N, 2)
  e <- matrix(rnorm(N, 1), N, 1)
  beta <- c(7, 3)

  y <- x%*% beta + 5 + e
  y <- y + sapply(inds[,1], function(company) company.fe[[company]])
  y <- y + sapply(inds[,2], function(person) person.fe[[person]])

  mget(c("y", "x", "e"), environment())
}

list2env(make.small.data(), environment())

test_that("Demean with fixed effects consisting of multiple disconnected components should be warned.", {
  expect_warning(solveFE(cbind(y, x), inds), warning.contain.multiple.components(), fixed = TRUE)
})

model <- suppressWarnings(solveFE(cbind(y, x), inds))

test_that("Predict with indicators from disconnected components should be warned", {
  N <- 3
  x <- matrix(rnorm(N * 2, 3), N, 2)
  inds <- rbind(c("facebook", "david"),
                c("facebook", "bob"),
                c("google", "alice"))

  warnings <- c()
  suppressWarnings(withCallingHandlers(predictFE(model, x, inds),
    warning = function(e) { warnings <<- c(warnings, conditionMessage(e)) }))

  # note that `factor()` will sort companies by their names
  error1 <- list(row = 1, group.x = 1, group.y = 2, value.x = 2, value.y = 4)
  error2 <- list(row = 3, group.x = 1, group.y = 2, value.x = 3, value.y = 1)

  expect_equal(warnings[1], warning.cross.component(error1, inds, model$group.levels))
  expect_equal(warnings[2], warning.cross.component(error2, inds, model$group.levels))

  expect_equal(length(warnings), 2)
})

test_that("Predict with indicators from disconnected components should yield NA", {
  N <- 3
  x <- matrix(rnorm(N * 2, 3), N, 2)
  inds <- rbind(c("facebook", "david"),
                c("facebook", "bob"),
                c("google", "alice"))

  predicted <- suppressWarnings(predictFE(model, x, inds))
  expect_equal(predicted[1], as.numeric(NA))
  expect_equal(predicted[3], as.numeric(NA))
})
