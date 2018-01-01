SEED  <- 19260817
N     <- 2000
LEVEL <- 50

make.small.data <- function () {
  set.seed(SEED)

  x <- matrix(rnorm(N * 5, 3), N, 5)
  e <- matrix(rnorm(N, 1), N, 1)

  inds <- matrix(sample(LEVEL, N * 3, replace = TRUE), N, 3)
  effs <- matrix(runif(LEVEL * 3), LEVEL, 3)

  with.effects <- function(j) sapply(1 : N, function(i) effs[inds[i, j], j])

  y <- 7 * x[,1] + 3 * x[,2] + 2 * x[,3] + 5 * x[,4] + 8 * x[,5] + 5 + e +
       with.effects(1) + with.effects(2) + with.effects(3)

  mget(c("y", "x", "e", "inds", "effs"), environment())
}

change.group.indicators <- function (inds) {
  prefixes <- c("alpha", "beta", "gamma")
  inds <- sapply(1 : 3, function(i) {
    rename <- function(j) paste(prefixes[[i]], inds[j, i], sep = ".")
    sapply(1 : N, rename)
  })
}
