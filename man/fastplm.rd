\name{fastplm}
\alias{fastplm}
\title{Solve Fixed Effect Model}

\description{\code{fastplm} solves a fixed effects model or a 2SLS regression with fixed effects.}

\usage{fastplm(data = NULL, formula = NULL,
			   index = NULL, y = NULL, 
			   x = NULL, z = NULL,
			   ind = NULL, sfe = NULL,
			   cfe = NULL, PCA = TRUE,
			   sp = NULL, knots = NULL,
			   degree = 3, se = 1,
			   vce = "robust", cluster = NULL,
			   wild = TRUE, refinement = FALSE,
			   test_x = NULL, parallel = TRUE,
			   nboots = 200, seed = NULL,
			   core.num = 1, drop.singletons = TRUE,
			   bootcluster = NULL, iv.test = TRUE,
			   first = FALSE, orthog = NULL,
			   endog = NULL)}

\arguments{

  \item{data}{An object of class "dataframe" or "matrix". If both \code{formula} 
   and \code{y} are omitted, the first column of \code{data} will be regarded 
   as the vector of outcome variable and other columns will be regarded as 
   regressors. If both \code{data} and \code{y} present, \code{data} is omitted.}
   
  \item{formula}{An object of class "formula": a symbolic description of
   the model to be fitted. If it has the form "Y~1", 
   the function merely estimates fixed effects.}
  
  \item{index}{A string vector specifying the indicators. Omissible if 
   \code{formula} is omitted.}

  \item{y}{The Y matrix. Omissible if \code{formula} is provided.}

  \item{x}{The X matrix. Omissible if \code{formula} is provided. 
   If \code{X} are not provided and \code{formula} is absent, 
   the model is still valid. In such case, 
   the function merely estimates fixed effects. }
   
  \item{z}{A string vector of the names of all exogenous variables 
  (both of excluded IVs and included IVs) if \code{formula} is present. 
   It can also be a matrix of all exogenous variables if 
   \code{formula} is absent. }

  \item{ind}{A matrix where each row corresponds to the row (observation) in the 
    X-Y dataset and each column represents an effect. The entry of row X of effect 
    Y represents the group in Y that X belongs to.
    Groups can be specified by either numbers or strings, 
	i.e. the input matrix can be of mode \code{numeric} or \code{character}. }

  \item{sfe}{A vector index of simple (additive) fixed effects. 
    Each can be specified either by effect name or position as in \code{index} or 
    by the position in the indicator matrix \code{ind}.
	If omitted, \code{sfe} will be the collection of all effects in 
    \code{index} or \code{inds}. }

  \item{cfe}{A list index of complex fixed effects. 

    Specifically, an complex fixed effect is a generalized fixed effect. 
	It consists a pair of two effects, (I, E), interacting with each other. 
	"I" stands for influence, whose level has an observed vector \code{weight}. 
	"E" stands for effect, whose level has an unobserved vector \code{coefficient} 
	to estimate. For each observation, 
	the effect of (I, E) is the dot product of the weight of the row's 
    level for "I" and the coefficient of the row's level for "E".  

    Each element is a vector whose length is 2. The 1st item of each element is the 
    index of effect and the 2nd item is the index of influence. For index, see 
    \code{sfe}.

    If omitted, complex fixed effects will not be estimated. 

  }

  \item{PCA}{A logical flag indicating whether to perform principal 
	components analysis for influence in complex fixed effects (see \code{cfe}).}

  \item{sp}{A character value or a numeric vector specifying the variable for 
    fitting a b-spline curve.}

  \item{knots}{A numeirc value specifying the knots point for \code{sp}. If 
    left blank, a polynomial curve will be fitted by default.}

  \item{degree}{A positive integer speficying the order of the spline curve. 
    Default value is \code{degree = 3} for a cubic curve.}

  \item{se}{A logical flag indicating whether uncertainty estimates of 
    covariates will be produced. For 2sls regressions, all diagnostic tests will
	not be conducted if \code{se = FALSE}.}

  \item{vce}{A character value indicating type for variance estimator. 
    Choose from: "standard" for standard ols standard errors, "robust" for 
    the Huber White robust standard errors (default value), "clustered" (or "cl") 
    for clustered standard errors, "jackknife" for jackknife standard errors, 
    and "bootstrap" (or "boot") for bootstrapped standard errors.}

  \item{cluster}{A character value of the clustered variable(s) in the data frame if 
    \code{formula} is provided or a matrix object of the clustered variable(s) for 
    clustered standard error or clustered bootstrap. 
	Two-way clustering is also supported. }

  \item{wild}{A logical flag specifies if wild bootstrap will be performed to obtain 
    uncertainty estimates. Omissible if \code{se = FALSE}. 
	Default to TRUE for computational efficiency.}

  \item{refinement}{A logical flag specifies if bootstrap refinement will be 
    performed to obtain uncertainty estimates. Omissible if \code{se = FALSE}.
	Default to FALSE.}  

  \item{test_x}{A character specifies the variable of interest for wild 
    bootstrap refinement. Omissible if \code{refinement = FALSE} or 
    \code{wild = FALSE}. When \code{wild = TRUE}, if \code{test_x = NULL}, 
	an unrestricted wild bootstrap will be conducted, otherwise, 
	a restricted wild bootstrap will be performed. Default to NULL.}  

  \item{parallel}{A logical flag indicating whether to perform parallel 
	computing for the bootstrap procedure. Default to TRUE.}

  \item{nboots}{An integer specifying the number of bootstrap
    runs. Omissible if \code{se = FALSE}.}

  \item{seed}{An integer that sets the seed in random number generation. 
    Omissible if \code{se = FALSE}.}

  \item{core.num}{The number of cores that will be used for computation. 
	Default is one. Larger number of cores can improve the speed moderately.}

  \item{drop.singletons}{A logical flag indicating whether to drop 
    singletons(groups with only one observation). Default to TRUE.
	The program can report the number of singletons it detects and return a vector
	\code{drop.index} noting the indices of dropped observations.}
	
  \item{bootcluster}{A character value specifying the level of bootstrap in 
	two-way clustered bootstrap refinement. 
	Must be one of the names of \code{cluster}.}

  \item{iv.test}{A logical flag indicating whether to conduct diagnostic tests 
	for 2SLS regressions, including weak IV tests, exogeneity tests, 
	and endogeneity tests. Default to TRUE. Omissible if \code{se=FALSE} or 
	\code{vce="bootsrap" or "jackknife"}.}
	
  \item{first}{A logical flag indicating whether to compute partial F statistics,
	Sanderson-Windmeijer (SW) F statistics and Angrist and Pischke (AP) F statistics
	for first stage regression for 2SLS estimation. Default to FALSE.}
	
  \item{orthog}{A character value(s) specifying the exogenous variables whose exogeneity
   is to be tested. A C-statistic will be calculated in this case. Default to NULL.}
  
  \item{endog}{A character value(s) specifying the endogenous variables whose endogeneity
   is to be tested. A C-statistic will be calculated in this case. Default to NULL.}

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
  
  \item{est.coefficients}{The estimated coefficients, standard errors and confidence intervals.}

  \item{vcov}{The estimated variance-covariance matrix of estimated coefficients.}
  
  \item{R2}{Multiple R-squared(original model).}
  
  \item{Adj_R2}{Adjusted multiple R-squared(original model).}
  
  \item{projR2}{Multiple R-squared(demeaned model).}
  
  \item{projAdj_R2}{Adjusted multiple R-squared(demeaned model).}
  
  \item{RMSE}{Root Mean Square Error.}
  
  \item{proj_F_statistic}{F-statistic.}
  
  \item{variance.type}{Type of standard errors and vcov matrices.}
  
  \item{num.singletons}{Number of singletons(groups with only one observation) the program detected.}
  
  \item{drop.index}{The original index of dropped singletons.}

  \item{inds}{As in the input.}

  \item{fe}{An object of "fixed.effects" created for intermediate computing.}

  \item{df.residual}{Degree of freedom.}
  
  \item{dof}{A list that stores the number of categories and the number of 
	redundant categories in each level of fixed effects.}

  \item{refinement}{A list that restores results from wild/pairwise bootsrap refinement.}

  \item{tests}{Only for 2SLS regressions, stores the results of diagnostic tests.}
  
  \item{ols.model}{Only for 2SLS regressions, stores the results of the simple 
	ols regression (only regress y on x).}
	
  \item{reduce.model}{Only for 2SLS regressions, stores the results of the reduced-form
	ols regression (regress y on all exogenous variables z).}
	
  \item{names}{Only for 2SLS regressions, stores the names of excluded instruments and 
	included instruments.}

}

\examples{
data(fastplm)
out <- fastplm(y~x1+x2, index = c("country", "industry", "year"),
               data = simdata1, 
               se = TRUE)
}


