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

  beta <- c(7, 3, 2, 5, 8)

  y <- x %*% beta + 5 + e + rowSums(sapply(1 : 3, with.effects))

  mget(c("y", "x", "e", "inds", "effs"), environment())
}

change.group.indicators <- function (inds) {
  prefixes <- c("alpha", "beta", "gamma")
  inds <- sapply(1 : 3, function(i) {
    rename <- function(j) paste(prefixes[[i]], inds[j, i], sep = ".")
    sapply(1 : N, rename)
  })
  colnames(inds) <- prefixes
  return(inds)
}
