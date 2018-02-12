#ifndef FASTPLM_FIXED_EFFECTS_H
#define FASTPLM_FIXED_EFFECTS_H

#include "Common.h"
#include "ComponentAnalysis.h"

struct SimpleFixedEffect {
    std::size_t size;
    const arma::uvec levelSizes;
    const arma::subview_col<arma::uword> indicators;

    SimpleFixedEffect(std::size_t size, const arma::uvec levelSizes, const arma::subview_col<arma::uword> indicators):size(size), levelSizes(levelSizes), indicators(indicators) {}

    arma::vec demean(arma::subview_col<double> data) const;
};

struct FixedEffects {
public:
    std::size_t size;
    arma::uvec groupSizes;
    arma::umat indicators;
    std::vector<const SimpleFixedEffect> simpleEffects;

    optional<ComponentTables> componentTables;
    std::vector<CrossComponentError> checkComponents(const arma::mat& indicators) const;

    std::vector<arma::mat> demean(arma::mat& data) const;
    static std::unique_ptr<const FixedEffects> create(const arma::uvec& groupSizes, const arma::mat& indicators);
};

#endif
