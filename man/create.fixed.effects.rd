\name{create.fixed.effects}
\title{Create Fixed Effect System}

\description{\code{create.fixed.effects} creates a fixed effect system for fitting.}

\usage{create.indicators(inds, sfes = NULL, cfes = NULL)}

\arguments{
  \item{inds}{An object of class \code{indicators} created by \code{create.indicators}.}

  \item{sfes}{A vector of simple (i.e. non-complex) fixed effects. Each can be specified either by effect name as in \code{inds$effect.names} or by the index of the effect in \code{inds}.

  If omitted, \code{sfes} will be the collection of all effects in \code{inds}.}

  \item{cfes}{A list of complex fixed effects. Each should be an object
    of class \code{complex.effect} created by \code{create.complex.effect}.

  If omitted, \code{cfes} will be the empty list.}
}

\value{An object of class \code{fixed.effects}.

Two properties to note are \code{sfe.names} and \code{cfe.names}. They represent the names of those effects. Names for simple fixed effects will be taken from \code{inds}, whereas names for complex fixed effects will be of the form \code{"effect.1:effect.2"}, where \code{"effect.1"} is the effect
and \code{"effect.2"} is the influence.}

\note{It is possible that the given \code{inds} may not be connected, which will be detected in \code{create.fixed.effects}. This is not an issue for fitting the model but could lead to meaningless prediction.

By connectedness, we mean the graph formed by \code{inds} is not connected. Consider the case where Alice and Bob have been employees of both Apple and Microsoft, and Cathy and David employees of both Google and Facebook. As a graph, there is no edge between the group Alice, Bob, Apple, and Microsoft, and the group Cathy, David, Google, and Facebook. If we add an arbitrary constant to the estimated effects for Alice and Bob and subtract that constant from the estimated effects for Apple and Microsoft, we will get a different yet equally valid fitted solution. In such case, however, cross-component predicition (e.g. Alice at Google) is meaningless.}

\author{Minsheng Liu}
