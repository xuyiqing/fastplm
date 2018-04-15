#ifndef FASTPLM_FIXED_EFFECTS_H
#define FASTPLM_FIXED_EFFECTS_H

#include "FixedEffect.h"
#include "ComponentAnalysis.h"

struct FixedEffects {
public:
    std::size_t size;
    arma::uvec levelCounts;
    arma::umat indicators;
    std::vector<const SimpleFixedEffect> simpleEffects;

    optional<ComponentTables> componentTables;
    std::vector<CrossComponentError> checkComponents(const arma::mat& indicators) const;

    std::vector<arma::mat> demean(arma::mat& data) const;
    static std::unique_ptr<const FixedEffects> create(const arma::uvec& levelCounts, const arma::mat& indicators);
};

#endif
