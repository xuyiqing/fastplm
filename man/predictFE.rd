\name{predictFE}
\alias{predictEF}
\title{Predict with a linear model with multi-way fixed effects fast}
\usage{predictFE(model, newX, FEValues = NULL, grandMean = 0)}

\description{Once a model is estimated by \code{solveFE}, \code{predictFE} can be used to predict the outcome given new data (including the grand mean, new X and new FE indicators).}

\arguments{
  \item{model}{The result of a call to \code{solveFE}.}
  
  \item{newX}{The new \eqn{X}.}

  \item{FEValues}{New group indicators. Note that if the argument is supplied, the result given in \code{model} must have its fixed effects estimated (see more at the documentation for \code{solveFE}.
  
  Warning. Currently we do not check the validity of supplied group indicators. Therefore, it is the caller's responsibility to make sure that all group indicators are contained in the original data passed to \code{solveFE}. This will change in a future release.}
  
  \item{grandMean}{An optional grandMean. It satisifies \code{predictFE(m, x, fe, mean) = mean + predictFE(m, x, fe, 0)}.}
}

\value{The outcome vector.}

\author{Minsheng Liu}

\examples{
set.seed(57)

n <- 10000
nlvl <- 40

x1 <- rnorm(n, 3)
x2 <- rnorm(n, 3)
x3 <- rnorm(n, 3)
x4 <- rnorm(n, 3)
e <- rnorm(n, 1) # error

## generate 3 group indicators
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

## execute with two cores
coreNum <- 2
result <- solveFE(as.matrix(d[,c("y","x1","x2","x3","x4")]), as.matrix(d[,c("gp1","gp2","gp3")]), coreNum)

outcome <- predictFE(result, as.matrix(d[,c("x1","x2","x3","x4")]), as.matrix(d[,c("gp1","gp2","gp3")]))
}
