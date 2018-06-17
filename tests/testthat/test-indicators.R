context("Indicators Related Functionality")

source("plain-data.R")

test_that("Creating indicators should be independent of concrete values of indicators.", {
  # NOTE: it is required that factor(raw.inds) and factor(new.inds)
  # gives the same result, i.e. the two indicator matrices are isomorphic
  # up to order-preserving renaming.

  inds1 <- create.indicators(raw.inds)
  inds2 <- create.indicators(new.inds)

  expect_equal(class(inds1), "indicators")

  expect_equal(inds1$level.sizes, inds2$level.sizes)
  expect_equal(inds1$inds, inds2$inds)

  expect_equal(inds1$effect.names, c("effect.1", "effect.2", "effect.3"))
  expect_equal(inds2$effect.names, colnames(new.inds))
})

test_that("Creating sub-indicaotrs should handle unknown indicators.", {
  row  <- sample(N, 1)
  col  <- sample(3, 1)
  inds <- create.indicators(new.inds)
  fe   <- create.fixed.effects(inds)

  our.inds <- new.inds
  our.inds[row, col] <- "unknown"

  inds1 <- expect_warning(create.subindicators(our.inds, inds = inds, fe = fe),
    WARN.INDS.unknown.indicator(row, col, our.inds),
    fixed = TRUE)

  inds2 <- create.subindicators(new.inds, inds = inds, fe = fe)

  actual          <- inds1$inds
  expected        <- inds1$inds
  expected[row, ] <- NA
  expect_identical(actual, expected)
})

test_that("Creating indicators/sub-indicaotrs should catch potential disconnected components", {
  company.fe <- as.list(rnorm(4, 10))
  names(company.fe) <- c("facebook", "microsoft", "amazon", "google")

  person.fe <- as.list(rnorm(4, 10))
  names(person.fe) <- c("alice", "bob", "cathy", "david")

  raw.inds <- rbind(c("facebook", "alice"),
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
    y <- y + sapply(raw.inds[,1], function(company) company.fe[[company]])
    y <- y + sapply(raw.inds[,2], function(person) person.fe[[person]])

    mget(c("y", "x", "e"), environment())
  }

  list2env(make.small.data(), environment())

  inds <- create.indicators(raw.inds)
  fe   <- expect_warning(create.fixed.effects(inds),
    WARN.FE.multiple.components(),
    fixed = TRUE)

  raw.sub.inds <- rbind(c("facebook", "david"),
                        c("facebook", "bob"),
                        c("google", "alice"))
  warnings <- c()
  sub.inds <- suppressWarnings(withCallingHandlers(create.subindicators(raw.sub.inds, inds = inds, fe = fe),
    warning = function(e) { warnings <<- c(warnings, conditionMessage(e)) }))

  # note that `factor()` will sort companies by their names
  err1 <- list(row = 1, group.x = 1, group.y = 2, value.x = 2, value.y = 4)
  err2 <- list(row = 3, group.x = 1, group.y = 2, value.x = 3, value.y = 1)

  expect_equal(warnings[1], WARN.INDS.cross.component(err1, raw.sub.inds, inds$levels))
  expect_equal(warnings[2], WARN.INDS.cross.component(err2, raw.sub.inds, inds$levels))

  expect_equal(length(warnings), 2)
})
