#ifndef FASTPLM_MAP_H
#define FASTPLM_MAP_H

#include "Common.h"
#include "FixedEffect.h"
#include "PlainModel.h"

struct FastFESolver: public PlainModel {
public:
    enum DemeanTransform {
        kHalperin, kSymmetricHalperin, kCimmino
    };
    
    FastFESolver(const arma::mat& X, const arma::colvec& Y,
                 const std::vector<FixedEffect>& fixedEffects,
                 DemeanTransform transform = kHalperin,
                 bool doesComputeFixedEffects = false);
    void compute();
    
private:
    DemeanTransform transform;
    
    const bool doesComputeFixedEffects;
    
    arma::colvec Y_;
    arma::mat X_;
    
    void demean(arma::mat& data, const FixedEffect& effect);
    void halperin(arma::mat& data);
    void symmetricHalperin(arma::mat& data);
    
    void demean();
    
    void estimateParams();
    void estimateFixedEffects();
};

#endif
