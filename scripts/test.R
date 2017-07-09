
install.packages("rbenchmark")
install.packages("lfe")
install.packages("Rcpp")

setwd("~/Github/fastplm")

## Simulating FE sample

set.seed(123)
n <- 1000000
nlvl <- 20
x1 <- rnorm(n, 3)
x2 <- rnorm(n, 3)
x3 <- rnorm(n, 3)
x4 <- rnorm(n, 3)
x5 <- rnorm(n, 3)
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
y <- 5 + 7 * x1 + 3 * x2 + 2 * x3 + 5 * x4 + 8 * x5 + gp.eff[,1] + gp.eff[,2] + gp.eff[,3] + e 

## put in a data frame
d <- cbind.data.frame("id"=1:n, gp, y, x1, x2, x3, x4, x5, gp.eff, e)

## speedtest
library(lfe)
library(rbenchmark)
benchmark(
    outFELM    <- felm(y ~ x1 + x2 + x3 + x4 + x5 | gp1 + gp2 + gp3, data = d),
    outFastPLM <- fastplm(as.matrix(d[,c("y","x1","x2","x3","x4","x5")]), as.matrix(d[,c("gp1","gp2","gp3")])),
    order = NULL,
    replications = 1
)

## results
## replications elapsed relative user.self sys.self user.child sys.child
## 1 ols      100   0.789    1.672     0.783    0.005          0         0
## 2 lfe      100   6.840   14.492     1.613    0.048          0         0
## 3 ours     100   0.472    1.000     0.472    0.001          0         0

