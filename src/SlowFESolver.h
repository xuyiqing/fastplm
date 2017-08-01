#ifndef FASTPLM_SLOW_FESOLVER_H
#define FASTPLM_SLOW_FESOLVER_H

#include "Common.h"
#include "FixedEffect.h"
#include "PlainModel.h"

struct SlowFESolver: public PlainModel {
private:
    arma::colvec nextParams(const arma::colvec&);
    arma::colvec nextFixedEffect(const FixedEffect&, const arma::colvec&);
    double nextIntercept(const arma::colvec&);
    arma::colvec computeResidual(const Result&);

public:
    SlowFESolver(const arma::mat& X, const arma::colvec& Y,
                    const std::vector<FixedEffect>& fixedEffects,
                    bool computeXtXInverse = true);
    
    void compute();
};

#endif
