#ifndef FASTPLM_LINEAR_MODEL_H
#define FASTPLM_LINEAR_MODEL_H

#include "Common.h"

struct LinearModel {
public:
    arma::mat X;
    arma::vec Y;

    bool isLinearDependent;
    arma::uvec dependents, independents;
    arma::vec beta;

    static const LinearModel solve(const arma::mat& X, const arma::vec& Y);
    static const LinearModel solve(const arma::mat& data);

#ifndef BUILD_WITHOUT_R
    operator Rcpp::List() const {
        Rcpp::List _;
        _["x"] = X;
        _["y"] = Y;
        _["is.linear.dependent"] = isLinearDependent;
        _["dependents"] = dependents;
        _["independents"] = independents;
        _["beta"] = beta;
        return _;
    }
#endif
};

#endif
