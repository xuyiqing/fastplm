\name{solve.fe.model}
\title{Solve Fixed Effect Model}

\description{\code{solve.fe.model} solves a fixed effects model.}

\usage{solve.fe.model(inds, fe = NULL, data = NULL, x = NULL, y = NULL, core.num = 1)}

\arguments{
  \item{inds}{An object of class \code{indicators} created by \code{create.indicators}.}

  \item{fe}{An object of class \code{fixed.effects} created by \code{create.fixed.effects}.

  If omitted, \code{fe = create.fixed.effects(inds)}. The default case is all effects in \code{inds} are simple fixed effects, and there is no complex fixed effect.}

  \item{data}{The X-Y dataset. The first column is the Y vector.
    The rest is the X matrix.

    If omitted, \code{data = cbind(y, x)}. If both present, \code{data} is used.}

  \item{x}{The X matrix. Omissible if \code{data} is provided.

  If both \code{data} and \code{X} are not provided, the model is still valid. In such case, the function merely estimates fixed effects.}

  \item{y}{The Y vector. Omissible if \code{data} is provided.}

  \item{core.num}{The number of cores that will be used for computation. Default is one.

  Do not use more than the number of your physical cores.}
}

\value{An object of class \code{fe.model} with at least the following properties:

  demeaned: A list represents the demeaned linear model. It has, among other properties, \code{x} and \code{y}, which are the demean result of the input \code{x} and \code{y}.

  coefficients: The coefficients for the demeaned linear model.

  sfe.coefs: A list of estimated simple fixed effects. Each is a column vector, in which row names correspond to names of levels in \code{inds}. Each, as a list item, has the same name as its effect name in \code{inds}. Each vector is centered, i.e. set to have zero mean.

  cfe.coefs: A list of estimated complex fixed effects. Each is a matrix in which each row represents a level. The naming convention is the same as in \code{sfe.coefs}. There is no centering. Expect an arbitrary intercept.

  fitted.values: The fitted values of the fixed effect model.

  residuals: The residuals of the fixed effect model, which is numerically the same as the residual of the demeaned linear model.

  intercept: The intercept of the fixed effect model. Arbitrary if complex fixed effects are present.

  inds: As in the input.

  fe: As in the input (or created, as documented above).
}

\examples{
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

###########################################

inds  <- create.indicators(raw.inds)
model <- solve.fe.model(inds, y = y, x = x)
}

\author{Minsheng Liu}
