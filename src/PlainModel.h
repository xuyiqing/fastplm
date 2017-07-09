#ifndef FASTPLM_PLAIN_MODEL_H
#define FASTPLM_PLAIN_MODEL_H

#include "Common.h"
#include "FixedEffect.h"

struct PlainModel {
public:
    arma::uword paramCount;
    arma::uword observationCount;
    
    const arma::mat X;
    const arma::colvec Y;
    const std::vector<FixedEffect>& fixedEffects;
    
    const arma::mat XtXInverse;
    
    PlainModel(const arma::mat& X, const arma::colvec& Y,
               const std::vector<FixedEffect>& fixedEffects,
               bool computeXtXInverse = true);
    
    struct Result {
        arma::colvec params;
        std::vector<arma::colvec> effects;
        double intercept;
        
        double computeDiff(const Result& last);
    };
    Result result;
};

#endif
