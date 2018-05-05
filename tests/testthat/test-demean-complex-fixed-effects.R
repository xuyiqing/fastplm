context("Demean Complex Fixed Effects")

SEED   <- 19260817
N      <- 2000
LEVELS <- c(100, 20)

set.seed(SEED)

x <- matrix(rnorm(N * 5, 3), N, 5)
e <- matrix(rnorm(N, 1), N, 1)

inds <- sapply(1 : 2, function(i) sample(LEVELS[i], N, replace = TRUE))

x[, 1] <- x[, 1] + inds[, 2]

sfes <- lapply(1 : 2, function(i) matrix(rnorm(LEVELS[i], 3), LEVELS[i], 1))
cfes <- lapply(1 : 2, function(i) matrix(rnorm(LEVELS[i], 3), LEVELS[i], 1))

inf1 <- rnorm(LEVELS[1])
inf2 <- (1 : LEVELS[2]) * 20
infs <- list(inf1, inf2)

with.sfe <- function(eff) sapply(1 : N, function(row) {
  level <- inds[row, eff]
  (sfes[[eff]])[level]
})

with.cfe <- function(inf, eff) sapply(1 : N, function(row) {
  levelEff <- inds[row, eff]
  levelInf <- inds[row, inf]
  (cfes[[eff]])[levelEff] * (infs[[inf]])[levelInf]
})

mapped.inf <- function(inf) sapply(1 : N, function(row) {
  level <- inds[row, inf]
  (infs[[inf]])[level]
})

sum.sfe <- rowSums(sapply(1 : 2, with.sfe))
sum.cfe <- with.cfe(1, 2) + with.cfe(2, 1)

beta <- c(7, 3, 2, 5, 8)

y <- x %*% beta + 5 + e + sum.sfe + sum.cfe

data <- cbind(y, x)
indicators <- create.indicators(inds)
cfe1 <- create.complex.effect(indicators, 1, 2, t(as.matrix(inf2)))
cfe2 <- create.complex.effect(indicators, 2, 1, t(as.matrix(inf1)))
cpp  <- create.fixed.effects(indicators, cfes = list(cfe1, cfe2))

result.lm      <- lm(y ~ x + factor(inds[, 1]) * mapped.inf(2)
                           + factor(inds[, 2]) * mapped.inf(1))
result.fastplm <- SolveFixedEffects(data, cpp, 1)

test_that("Demean complex fixed effects.", {
  expect_equal(as.matrix(unname(result.lm$coefficients[2 : 6])),
               unname(result.fastplm$coefficients))
})
