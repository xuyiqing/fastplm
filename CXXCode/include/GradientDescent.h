#ifndef FASTPLM_GRADIENT_DESCENT_H
#define FASTPLM_GRADIENT_DESCENT_H

#include "Common.h"
#include "FixedEffect.h"

struct GradientDescent {
    arma::uword paramCount;
    arma::uword observationCount;
    
    const arma::mat X;
    const arma::colvec Y;
    const std::vector<FixedEffect>& fixedEffects;
    
    const arma::mat XtXInverse;
    
    GradientDescent(const arma::mat& X_, const arma::colvec& Y_, const std::vector<FixedEffect>& fixedEffects);
    
    struct Result {
        arma::colvec params;
        std::vector<arma::colvec> effects;
        double intercept;
        
        double computeDiff(const Result& last);
    };
    Result result;
    
    arma::colvec nextParams(const arma::colvec&);
    arma::colvec nextFixedEffect(const FixedEffect&, const arma::colvec&);
    double nextIntercept(const arma::colvec&);
    
    arma::colvec computeResidual(const Result&);
    void compute();
};

#endif
