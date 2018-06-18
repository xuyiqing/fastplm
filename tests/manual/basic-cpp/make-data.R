SEED  <- 19260817
N     <- 64000
LEVEL <- 64000

set.seed(SEED)

x <- matrix(rnorm(N * 5, 3), N, 5)
e <- matrix(rnorm(N, 1), N, 1)

inds <- matrix(sample(LEVEL, N * 3, replace = TRUE), N, 3)
effs <- matrix(runif(LEVEL * 3), LEVEL, 3)

with.effects <- function(j) sapply(1 : N, function(i) effs[inds[i, j], j])

beta <- c(7, 3, 2, 5, 8)

y <- x %*% beta + 5 + e + rowSums(sapply(1 : 3, with.effects))

data <- cbind(y, x)

write.table(data, file = "data.csv", sep=",", row.names = FALSE, col.names = FALSE)
write.table(inds, file = "inds.csv", sep=",", row.names = FALSE, col.names = FALSE)
