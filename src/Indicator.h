#ifndef FASTPLM_INDICATOR_H
#define FASTPLM_INDICATOR_H

#include "Common.h"

struct Indicator {
    const arma::uvec indicator;
    const std::size_t levelCount;
    const arma::vec levelSizes;
    
    Indicator(arma::uvec&& indicator, std::size_t levelCount, arma::vec&& levelSizes): indicator(std::move(indicator)), levelCount(levelCount), levelSizes(std::move(levelSizes)) {}
};

std::vector<Indicator> createIndicators(const arma::uvec& levelCounts, const arma::mat& indsR);

#endif
