#ifndef FASTPLM_FIXED_EFFECTS_H
#define FASTPLM_FIXED_EFFECTS_H

#include "FixedEffect.h"
#include "ComponentAnalysis.h"

struct FixedEffects {
public:
    std::vector<Indicator> indicators;
    std::vector<SimpleFixedEffect> simpleEffects;
    std::vector<ComplexFixedEffect> complexEffects;

    optional<ComponentTables> componentTables;
    std::vector<CrossComponentError> checkComponents(const arma::mat& indicators) const;

    using Deltas = std::pair<std::vector<arma::vec>, std::vector<arma::mat>>;

    std::vector<Deltas> demean(arma::mat& data) const;
    static std::unique_ptr<const FixedEffects> create(const arma::uvec& levelCounts, const arma::mat& indsR, const arma::uvec& simpleEffects, const arma::uvec& complexEffects, const arma::uvec& complexInfluences, const std::vector<arma::mat>& weights);
};

#endif
