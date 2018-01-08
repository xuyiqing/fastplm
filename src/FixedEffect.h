#ifndef FASTPLM_FIXED_EFFECT_H
#define FASTPLM_FIXED_EFFECT_H

#include <unordered_map>

#include "Common.h"

struct FixedEffect {
    arma::uword groupCount;
    std::vector<int> groupSizes;
    std::vector<int> column;

    static FixedEffect fromColumn(const arma::colvec& column);
    void demean(arma::mat& data, arma::mat& deltas) const;

    std::vector<double> computeMean(const double* ptr) const;
private:
    FixedEffect();
};

#endif
