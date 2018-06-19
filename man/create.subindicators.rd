\name{create.subindicators}
\title{Create Subindicators}

\description{\code{create.subindicators} creates an object of class \code{subindicators} with levels from a given indicators, checked by the latter's fixed effect system. The result will be used in \code{predict.fe.model}.}

\usage{create.subindicators(sub.inds, model = NULL, inds = NULL, fe = NULL)}

\arguments{
  \item{sub.inds}{A new matrix descible the new observations for prediction, in the same format as the input for \code{create.indicators}.}

  \item{model}{The model solved by \code{solve.fe.model}. Its \code{inds} and \code{fe} will be used.

  If omitted, \code{inds} and \code{fe} must be provided explicitly. If both are present, \code{model} will be used.}

  \item{inds}{See \code{model} above.}

  \item{fe}{See \code{model} above.}
}

\value{An object of class \code{subindicators}.}

\note{As noted in the documentation for \code{create.fixed.effects}, there is a potential issue of disconnected indicators. As detailed there, the issue does not effect estimation but prediction. If cross-component observation is present in \code{sub.inds}, it will be warned, and the relevant rows will be set to \code{NA} and not predicted.}

\author{Minsheng Liu}
