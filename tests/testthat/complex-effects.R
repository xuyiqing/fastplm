SEED   <- 19260817
N      <- 2000
LEVELS <- c(100, 20)

set.seed(SEED)

x <- matrix(rnorm(N * 5, 3), N, 5)
e <- matrix(rnorm(N, 1), N, 1)

raw.inds <- sapply(1 : 2, function(i) sample(LEVELS[i], N, replace = TRUE))

x[, 1] <- x[, 1] + raw.inds[, 2]

sfes <- lapply(1 : 2, function(i) matrix(rnorm(LEVELS[i], 3), LEVELS[i], 1))
cfes <- lapply(1 : 2, function(i) matrix(rnorm(LEVELS[i], 3), LEVELS[i], 1))

inf1 <- rnorm(LEVELS[1])
inf2 <- (1 : LEVELS[2]) * 20
infs <- list(inf1, inf2)

with.sfe <- function(eff) sapply(1 : N, function(row) {
  level <- raw.inds[row, eff]
  (sfes[[eff]])[level]
})

with.cfe <- function(inf, eff) sapply(1 : N, function(row) {
  level.eff <- raw.inds[row, eff]
  level.inf <- raw.inds[row, inf]
  (cfes[[eff]])[level.eff] * (infs[[inf]])[level.inf]
})

mapped.inf <- function(inf) sapply(1 : N, function(row) {
  level <- raw.inds[row, inf]
  (infs[[inf]])[level]
})

sum.sfe <- rowSums(sapply(1 : 2, with.sfe))
sum.cfe <- with.cfe(1, 2) + with.cfe(2, 1)

beta <- c(7, 3, 2, 5, 8)

y <- x %*% beta + 5 + e + sum.sfe + sum.cfe
