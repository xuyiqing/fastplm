#ifndef FASTPLM_FAST_FESOLVER_H
#define FASTPLM_FAST_FESOLVER_H

#include "Common.h"
#include "FixedEffect.h"
#include "LinearModel.h"

struct FixedEffectModel {
public:
    LinearModel demeaned;

    arma::vec coefficients;
    std::vector<arma::vec> feCoefs;

    arma::vec fittedValues, residuals;
    double intercept;

    static const FixedEffectModel solve(const arma::mat& data, const std::vector<FixedEffect>& fixedEffects);
    static const FixedEffectModel solve(const arma::mat& X, const arma::vec& Y, const std::vector<FixedEffect>& fixedEffects) {
        return solve(arma::join_horiz(Y, X), fixedEffects);
    }

#ifndef BUILD_WITHOUT_R
    operator Rcpp::List() const {
        Rcpp::List _;
        _["demeaned"] = static_cast<Rcpp::List>(demeaned);
        _["coefficients"] = demeaned.beta;

        Rcpp::List feCoefs_;
        for (const auto& coefs : feCoefs)
            feCoefs_.push_back(coefs);
        _["FEcoefs"] = feCoefs_;

        _["fitted.values"] = fittedValues;
        _["residuals"] = residuals;
        _["intercept"] = intercept;

        return _;
    }
#endif
};

#endif
