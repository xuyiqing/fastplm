SEED  <- 19260817

make.data <- function(N, LEVEL) {
  set.seed(SEED)

  x <- matrix(rnorm(N * 5, 3), N, 5)
  e <- matrix(rnorm(N, 1), N, 1)

  inds <- matrix(sample(LEVEL, N * 3, replace = TRUE), N, 3)
  effs <- matrix(runif(LEVEL * 3), LEVEL, 3)

  with.effects <- function(j) sapply(1 : N, function(i) effs[inds[i, j], j])

  beta <- c(7, 3, 2, 5, 8)

  y <- x %*% beta + 5 + e + rowSums(sapply(1 : 3, with.effects))

  list(x = x, y = y, inds = inds)
}

library(lfe)
library(rbenchmark)
library(fastplm)

run <- function(N) {
  data <- make.data(N, 100)
  x    <- data$x
  y    <- data$y
  inds <- data$inds

  benchmark(
    R1 <- felm(y ~ x | inds[, 1] + inds[, 2] + inds[, 3]),
    R2 <- solve.fe.model(create.indicators(inds), x = x, y = y),
    order = NULL,
    replications = 100
  )
}
