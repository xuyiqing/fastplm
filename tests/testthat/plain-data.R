SEED  <- 19260817
N     <- 2000
LEVEL <- 50

set.seed(SEED)

x <- matrix(rnorm(N * 5, 3), N, 5)
e <- matrix(rnorm(N, 1), N, 1)

raw.inds <- matrix(sample(LEVEL, N * 3, replace = TRUE), N, 3)
sfe.coefs <- matrix(runif(LEVEL * 3), LEVEL, 3)

with.effects <- function(j) sapply(1 : N,
  function(i) sfe.coefs[raw.inds[i, j], j])

beta <- c(7, 3, 2, 5, 8)
effs <- rowSums(sapply(1 : 3, with.effects))

y <- x %*% beta + 5 + e + effs

prefixes <- c("alpha", "beta", "gamma")
new.inds <- sapply(1 : 3, function(i) {
  rename <- function(j) sprintf("%s.%02d", prefixes[[i]], raw.inds[j, i])
  sapply(1 : N, rename)
})
colnames(new.inds) <- prefixes
