#ifndef FASTPLM_MAP_H
#define FASTPLM_MAP_H

#include "Common.h"
#include "FixedEffect.h"
#include "PlainModel.h"

struct FastFESolver: public PlainModel {
public:
    FastFESolver(const arma::mat& X, const arma::colvec& Y,
                 const std::vector<FixedEffect>& fixedEffects,
                 bool doesComputeFixedEffects = false);
    void compute();
    
private:
    const bool doesComputeFixedEffects;
    
    arma::colvec Y_;
    arma::mat X_;
    
    void demean();
    
    void estimateParams();
    void estimateFixedEffects();
};

#endif
