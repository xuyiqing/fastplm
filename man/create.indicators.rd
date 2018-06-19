\name{create.indicators}
\title{Create Indicators}

\description{\code{create.indicators} creates an object describing to which group each observation belongs to for each fixed effect. The object will be used in other functions.}

\usage{create.indicators(inds)}

\arguments{
  \item{inds}{A matrix where each row corresponds to the row (observation) in the X-Y dataset and each column represents an effect. The entry of row X of effect Y represents the group in Y that X belongs to.

  Groups can be specified by either numbers or strings, i.e. the input matrix can be of mode \code{numeric} or \code{character}.}
}

\value{An object of class \code{indicators} with at least the following properties:

  levels: Distinct groups in each effect. It is a list where each item, corresponding to each column of the input matrix (effect), is a vector of levels (distinct groups) in that effect.

  The order of levels in each list item is not guaranteed.

  level.sizes: Total number of levels for each effect.

  effect.names: The name of each effect. Column names for the input matrix will be used if present. Otherwise, they will be named as \code{"effect.1"}, \code{"effect.2"}, and so on.
}

\author{Minsheng Liu}
