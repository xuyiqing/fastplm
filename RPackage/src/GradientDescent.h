#ifndef FASTPLM_GRADIENT_DESCENT_H
#define FASTPLM_GRADIENT_DESCENT_H

#include "Common.h"
#include "FixedEffect.h"
#include "PlainModel.h"

struct GradientDescent: public PlainModel {
private:
    arma::colvec nextParams(const arma::colvec&);
    arma::colvec nextFixedEffect(const FixedEffect&, const arma::colvec&);
    double nextIntercept(const arma::colvec&);
    arma::colvec computeResidual(const Result&);

public:
    GradientDescent(const arma::mat& X, const arma::colvec& Y,
                    const std::vector<FixedEffect>& fixedEffects,
                    bool computeXtXInverse = true);
    
    void compute();
};

#endif
