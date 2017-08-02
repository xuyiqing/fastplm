## Simulating FE sample

library(lfe)
library(rbenchmark)
library(fastplm)

testSpeed <- function(n, nlvl, CORE = 1) {
  set.seed(57)

  x1 <- rnorm(n, 3)
  x2 <- rnorm(n, 3)
  x3 <- rnorm(n, 3)
  x4 <- rnorm(n, 3)
  e <- rnorm(n, 1) # error

  ## generate 3 group indicators, each of 20 levels
  gp <- matrix(sample(1:nlvl, n*3, replace = TRUE), n, 3)
  colnames(gp) <- c("gp1","gp2","gp3")

  ## generate group effect, stored in a nlvl*3 matrix
  gp.coef <- matrix(runif(nlvl * 3),nlvl,3)

  ## assign group eff to each observation
  gp.eff <- matrix(NA, n, 3)
  for (i in 1:n) {
    for (j in 1:3) {
      gp.eff[i,j] <- gp.coef[gp[i,j],j]
    }
  }
  colnames(gp.eff) <- c("gp_eff1","gp_eff2","gp_eff3")

  ## outcome
  y <- 5 + 7 * x1 + 3 * x2 + 2 * x3 + 5 * x4 + gp.eff[,1] + gp.eff[,2] + gp.eff[,3] + e

  ## put in a data frame
  d <- cbind.data.frame("id"=1:n, gp, y, x1, x2, x3, x4, gp.eff, e)

  benchmark(
    R1 <- felm(y ~ x1 + x2 + x3 + x4 | gp1 + gp2 + gp3, data = d),
    R2 <- solveFE(as.matrix(d[,c("y","x1","x2","x3","x4")]), as.matrix(d[,c("gp1","gp2","gp3")]), CORE),
    order = NULL,
    replications = 10
  )
}

testSpeed(8000, 40, 2)
testSpeed(16000, 40, 2)
testSpeed(32000, 40, 2)
testSpeed(64000, 40, 2)
testSpeed(128000, 40, 2)
testSpeed(256000, 40, 2)
testSpeed(512000, 40, 2)
testSpeed(1024000, 40, 2)
