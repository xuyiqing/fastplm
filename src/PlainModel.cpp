#include "PlainModel.h"

PlainModel::PlainModel(const arma::mat& X, const arma::colvec& Y,
                       const std::vector<FixedEffect>& fixedEffects,
                       bool computeXtXInverse):
    paramCount(X.n_cols), observationCount(X.n_rows), X(X), Y(Y), fixedEffects(fixedEffects),
    XtXInverse(computeXtXInverse ? arma::inv(X.t() * X) :
               arma::mat(arma::zeros(paramCount, paramCount))) {}
