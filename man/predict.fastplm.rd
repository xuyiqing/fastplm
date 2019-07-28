\name{predict.fastplm}
\alias{predict.fastplm}
\title{Predict Fixed Effect Model}

\description{Predicts new Y for a given X based on a solved model.}

\usage{\method{predict}{fastplm}(object, data = NULL, x = NULL, ind, ...)}

\arguments{
  \item{object}{The model solved by \code{fastplm}.}

  \item{data}{A dataframe object.}

  \item{x}{New design matrix.}

  \item{ind}{New raw indicator matrix.}
  
  \item{...}{Other arguments.}
}

\value{Predicted Y.}

%\author{}
