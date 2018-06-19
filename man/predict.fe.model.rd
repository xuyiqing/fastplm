\name{predict.fe.model}
\title{Predict Fixed Effect Model}

\description{\code{predict.fe.model} predicts new Y for a given X based on a solved model.}

\usage{predict.fe.model(model, x, inds)}

\arguments{
  \item{model}{The model solved by \code{solve.fe.model}.}

  \item{model}{New X.}

  \item{inds}{New indicators describing X created by \code{create.subindicators}.

  Raw indicator matrix is not accepted. Furthermore, it will be checked that the input \code{inds} is indeed created from \code{model$inds} to avoid potential mistakes.}
}

\value{Predicted Y.}

\author{Minsheng Liu}
