
install.packages("rbenchmark")
install.packages("lfe")
install.packages("Rcpp")

setwd("~/Github/fastplm")

## Simulating FE sample

set.seed(123)
n <- 1000
nlvl <- 20
x1 <- rnorm(n, 3)
x2 <- rnorm(n, 3)
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
y <- 5 + 1 * x1 + 3 * x2 + gp.eff[,1] + gp.eff[,2] + gp.eff[,3] + e 

## put in a data frame
d <- cbind.data.frame("id"=1:n, gp, y, x1, x2, gp.eff, e)
head(d)

## Regressions

### ols
out1 <- lm(y~ x1 + x2 + as.factor(gp1) + as.factor(gp2) + as.factor(gp3), d)
coef(out1)[1:3]

## lfe
library(lfe)
out2 <- felm(y ~ x1 + x2 | gp1 + gp2 + gp3, data = d)
coef(out2) 

## fastplm
#library(Rcpp)
#sourceCpp("./cpp/fastplm.cpp")

#out3<-fastplm(as.matrix(d[,c("y","x1","x2")]), as.matrix(d[,c("gp1","gp2","gp3")]))
#out3$coefficients

## speedtest
library(rbenchmark)
benchmark(
    out1 <- lm(y~ x1 + x2 + as.factor(gp1) + as.factor(gp2) + as.factor(gp3), d),
    out2 <- felm(y ~ x1 + x2 | gp1 + gp2 + gp3, data = d),
    #out3<-fastplm(as.matrix(d[,c("y","x1","x2")]), as.matrix(d[,c("gp1","gp2","gp3")])),
    order = NULL
)

## results
## replications elapsed relative user.self sys.self user.child sys.child
## 1 ols      100   0.789    1.672     0.783    0.005          0         0
## 2 lfe      100   6.840   14.492     1.613    0.048          0         0
## 3 ours     100   0.472    1.000     0.472    0.001          0         0

