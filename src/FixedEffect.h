#ifndef FIXED_EFFECT_H
#define FIXED_EFFECT_H

#include <unordered_map>

#include "Common.h"

struct FixedEffect {
    arma::mat indicators;
    arma::uword groupCount;
    arma::uword observationCount;
    std::vector<int> valuesOccurences;
    
    typedef std::unordered_map<double, int> Indices;
    Indices indices;
    
    static FixedEffect fromColumn(const arma::colvec& column);
    void demean(arma::mat& data) const;
    
private:
    FixedEffect();
};

#endif
