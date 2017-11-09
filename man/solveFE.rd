\name{solveFE}
\alias{solveFE}
\title{Fit a linear model with multi-way fixed effects fast}
\usage{solveFE(data, fixedEffects, coreNum = 1, estimateFE = FALSE)}

\description{
  \code{solveFE} solves a linear model with multi-way fixed effects. It works similar to \code{felm} from the \code{felm} package, but performs much faster. Mathematically, we estimate \eqn{\mathbf{\beta}} of the following equation:
  
  \deqn{\mathbf{Y} = \mathbf{\beta} \mathbf{X} + \sum_i \mathbf{D_i}\mathbf{F_i} + \mathbf{\epsilon}}
  
  We use method of alternating projections (MAP) to remove fixed effects and then estimated the coefficients using OLS. Optionally, we can use gradient descent to estimate the remaining fixed effects \eqn{\mathbf{F_i}}.
  
  The algorithm employed here is the same as the one used in \code{felm}, but we manage to improve constant performance. For relatively small dataset (tens of thousands of rows), our function is more than 10x faster. For sufficiently large dataset (millions of rows), our package is still 4x as \code{felm}.
}

\arguments{
  \item{data}{An augmented matrix with its first column being \eqn{\mathbf{Y}} and remaining columns being \eqn{\mathbf{X}}.}
  
  \item{fixedEffects}{A matrix records the fixed effects. Its i-th column and j-th row determines which group of the category i does the row j belong to. For instance, if we have three individuals \eqn{1, 2, 3} and time point \eqn{2001, 2002} and organize the data into long form in a time major way (putting rows with the same time point together), we would have a fixed effect matrix of the following form:
  
  \deqn{\begin{bmatrix} 1 & 2001 \\ 2 & 2001 \\ 3 & 2001 \\ 1 & 2002 \\ 2 & 2002 \\ 3 & 2002 \\\end{bmatrix}}}

  \item{coreNum}{The number of CPU cores you want to commit to the computing tasks. Yes, we support parallelism! Two things to note. First, \code{coreNum} should never exceed the number of physical cores you have. Second, the speed improvement would not be propotional to the number of cores you have, as one task taking a significant amount of time is OLS, which is not computed in parallel here.
  
  The reason why we should use physical cores rather than virtual ones enabled by Intel's so-called hyper-threading technology is that the speed of our program depends more on cache locality (namely putting data inside the cache) than the raw speed of CPUs. However, hyper-threading merely doubles the number of registers inside each core, so it does not help but actually slow down things.}
  
  \item{estimateFixedEffects}{A boolean value that tells the program whether we will estimate fixed effects or not using gradient descent. Note that the algorithm is notoriously slow and takes forever to converge, so do not use it unless necessary.}
}

\value{
  \item{coefficients}{The estimated \eqn{\mathbf{\beta}} as a vector.}
  \item{intercept}{The intercept as a real number.}
  \item{fitted.values}{The fitted values, namely the result of the estimated linear function applied to the given \eqn{X}.}
  \item{residuals}{The residuals, namely the given \eqn{Y} minus the fitted values.}
  \item{FEcoefs}{A list of estimated fixed effects, named by column names of the input \code{fixedEffects}. If some \eqn{N}-th column is not named, the corresponding effect will be named as \code{effectN}. Only present if \code{estimateFixedEffects} is set \code{TRUE}.}
}

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
}
