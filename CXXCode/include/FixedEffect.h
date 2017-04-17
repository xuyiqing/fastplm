#ifndef FIXED_EFFECT_H
#define FIXED_EFFECT_H

#include <unordered_map>

#include "Common.h"

struct FixedEffect {
    arma::mat indicator;
    arma::uword groupCount;
    arma::uword observationCount;
    std::vector<int> valuesOccurences;
    
    typedef std::unordered_map<double, int> Indices;
    Indices indices;
    
    static FixedEffect fromColumn(const arma::colvec& column);
    
private:
    FixedEffect();
};

#endif
