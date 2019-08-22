\name{fastplm}
\alias{fastplm}
\title{Solve Fixed Effect Model}

\description{\code{fastplm} solves a fixed effects model.}

\usage{fastplm(formula = NULL, data = NULL, index = NULL, 
               y = NULL, x = NULL, ind = NULL, 
               sfe = NULL, cfe = NULL, 
               PCA = TRUE, se = 1, vce = "robust", 
               cluster = NULL, wild = FALSE,
               refinement = FALSE, test_x = NULL, parallel = FALSE, 
               nboots = 200, seed = NULL, core.num = 1)}

\arguments{

  \item{formula}{An object of class "formula": a symbolic description of
   the model to be fitted.}
  
  \item{data}{An object of class "dataframe" or "matrix". If both \code{formula} 
   and \code{y} are omitted, the first column of \code{data} will be regarded 
   as the vector of outcome variable. If both \code{data} and \code{y} present, 
   \code{data} is used.}

  \item{index}{A string vector specifying the indicators. Omissible if 
   \code{formula} is omitted.}

  \item{y}{The Y vector. Omissible if \code{data} is provided.}

  \item{x}{The X matrix. Omissible if \code{data} is provided. 
   If both \code{data} and \code{X} are not provided, the model is still valid. In 
   such case, the function merely estimates fixed effects. }

  \item{ind}{A matrix where each row corresponds to the row (observation) in the 
    X-Y dataset and each column represents an effect. The entry of row X of effect 
    Y represents the group in Y that X belongs to.

    Groups can be specified by either numbers or strings, i.e. the input matrix can 
    be of mode \code{numeric} or \code{character}. }

  \item{sfe}{A vector index of simple (i.e. non-complex) fixed effects. 
    Each can be specified either by effect name or position as in \code{index} or 
    by the position in the indicator matrix \code{ind}.

    If omitted, \code{sfe} will be the collection of all effects in 
    \code{index} or \code{inds}. }

  \item{cfe}{A list index of complex fixed effects. 

    Specifically, an complex fixed effect is a generalized fixed effect. It consists a 
    pair of two effects, (I, E), interacting with each other. "I" stands for influence, 
    whose level has an observed vector \code{weight}. "E" stands for effect, whose 
    level has an unobserved vector \code{coefficient} to estimate. For each 
    observation, the effect of I, E) is the dot product of the weight of the row's 
    level for "I" and the coefficient of the row's level for "E".  

    Each element is a vector whose length is 2. The 1st item of each element is the 
    index of effect and the 2nd item is the index of influence. For index, see 
    \code{sfe}.

    If omitted, complex fixed effects will not be estimated. 

  }

  \item{PCA}{A logical flag indicating whether to perform principal components analysis 
    for influence in complex fixed effects (see \code{cfe}).}

  \item{se}{A logical flag indicating whether uncertainty estimates of 
    covariates will be produced.}

  \item{vce}{A character value indicating type for variance estimator. 
    Choose from: "standard" for standard ols standard errors, "robust" for 
    the Huber White robust standard errors (default value), "clustered" (or "cl") 
    for clustered standard errors, "jackknife" for jackknife standard errors, 
    and "bootstrap" (or "boot") for bootstrapped standard errors.}

  \item{cluster}{A character value of the clustered variable(s) in the data frame if 
    \code{formula} is provided or a matrix object of the clustered variable(s) for 
    robust standard error. Two-way clustering is also supported. }

  \item{wild}{A logical flag specifies if wild bootstrap will be performed to obtain 
    uncertainty estimates. Omissible if \code{se = FALSE}.}

  \item{refinement}{A logical flag specifies if clutser bootstrap refinement will be 
    performed to obtain uncertainty estimates. Omissible if \code{se = FALSE}.}  

  \item{test_x}{A character specifies the variable of interest for wild cluster 
    bootstrap refinement. Omissible if \code{refinement = FALSE} or 
    \code{wild = FALSE}.}  

  \item{parallel}{A logical flag indicating whether to perform parallel computing for 
    the bootstrap procedure.}

  \item{nboots}{An integer specifying the number of bootstrap
    runs. Omissible if \code{se = FALSE}.}

  \item{seed}{An integer that sets the seed in random number generation. 
    Omissible if \code{se = FALSE}.}

  \item{core.num}{The number of cores that will be used for computation. Default is one.
    Do not use more than the number of your physical cores.}
}

\value{An object of class \code{fastplm} with at least the following properties:

  \item{demeaned}{A list represents the demeaned linear model. It has, among other 
    properties, \code{x} and \code{y}, which are the demean result of the input 
    \code{x} and \code{y}.}

  \item{coefficients}{The coefficients for the demeaned linear model.}

  \item{sfe.coefs}{A list of estimated simple fixed effects. Each is a column 
    vector, in which row names correspond to names of levels in \code{inds}. Each, as 
    a list item, has the same name as its effect name in \code{inds}. Each vector is 
    centered, i.e. set to have zero mean.}

  \item{cfe.coefs}{A list of estimated complex fixed effects. Each is a matrix in 
    which each row represents a level. The naming convention is the same as in 
    \code{sfe.coefs}. There is no centering. Expect an arbitrary intercept.}

  \item{fitted.values}{The fitted values of the fixed effect model.}

  \item{residuals}{The residuals of the fixed effect model, which is numerically the 
    same as the residual of the demeaned linear model.}

  \item{intercept}{The intercept of the fixed effect model. Arbitrary if complex fixed 
    effects are present.}

  \item{inds}{As in the input.}

  \item{fe}{An object of "fixed.effects" created for intermediate computing.}

  \item{refinement}{A list that restores results from cluster bootsrap refinement.}

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
model <- fastplm(y = y, x = x, ind = raw.inds, se = 0)
}

%\author{}
