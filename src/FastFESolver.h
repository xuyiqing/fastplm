#ifndef FASTPLM_FAST_FESOLVER_H
#define FASTPLM_FAST_FESOLVER_H

#include "Common.h"
#include "FixedEffect.h"
#include "PlainModel.h"

struct FastFESolver: public PlainModel {
public:
    FastFESolver(const arma::mat& X, const arma::colvec& Y,
                 const std::vector<FixedEffect>& fixedEffects);
    void compute();

    std::vector<arma::mat> deltas;

private:
    arma::colvec Y_;
    arma::mat X_;

    void demean();
};

struct Result {
    arma::colvec beta;
    arma::uvec dependents;
    arma::uvec independents;
    std::vector<arma::colvec> fixedEffects;
    arma::colvec residuals;
    arma::colvec fittedValues;
    double intercept;
};

Result estimate(const arma::mat& initData, const arma::mat& projData, const std::vector<arma::mat>& deltas);

#endif
