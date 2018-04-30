#include "Indicator.h"

std::vector<Indicator> createIndicators(const arma::uvec& levelCounts, const arma::mat& indsR) {
    std::vector<Indicator> indsC;
    
    for (auto i = 0u; i < indsR.n_cols; i ++) {
        const auto levelCount = levelCounts(i);
        arma::uvec ind(indsR.n_rows);
        arma::vec levelSizes = arma::zeros(levelCount);
        
        for (auto j = 0u; j < indsR.n_rows; j ++) {
            const auto x = static_cast<std::size_t>(indsR(j, i)) - 1;
            ind(j) = x;
            levelSizes[x] += 1.0;
        }
        
        indsC.emplace_back(std::move(ind), levelCount, std::move(levelSizes));
    }

    return indsC;
}
