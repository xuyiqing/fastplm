\name{create.complex.effect}
\title{Create Complex Effect}

\description{\code{create.complex.effect} creates a complex effect (an fixed effect involving interactions between two "simple" fixed effects). The object will be used to create a fixed effect system for fitting.

Specifically, an complex fixed effect is a generalized fixed effect. It consists a pair of two effects, (I, E), interacting with each other. "I" stands for influence, whose level has an observed vector \code{weight}. "E" stands for effect, whose level has an unobserved vector \code{coefficient} to estimate. FOr each observation, the effect of I, E) is the dot product of the weight of the row's level for "I" and the coefficient of the row's level for "E".}

\usage{create.indicators(inds, eff, inf, weight)}

\arguments{
  \item{inds}{An object of class \code{indicators} created by \code{create.indicators}.}

  \item{eff}{The effect of the complex fixed effect. It can be specified either by effect name as in \code{inds$effect.names} or by the index of the effect in \code{inds}.}

  \item{inf}{The influence of the complex fixed effect. It can be specified either by effect name as in \code{inds$effect.names} or by the index of the effect in \code{inds}.}

  \item{weight}{A matrix of weights associated with the influence. Each column represents the weight of the corresponding level. Weight should be given in the same order as the order of levels for influence in \code{inds$levels}.

  Only matrix is accepted even if the weight associated with each level is one-dimensional (i.e. a scalar). Use \code{as.matrix} when necessary.}
}

\value{An object of class \code{complex.effect}.}

\author{Minsheng Liu}
