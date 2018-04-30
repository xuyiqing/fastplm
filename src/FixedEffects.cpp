#include <queue>

#include "CrushQueue.h"
#include "FixedEffects.h"

arma::umat convertNumberingRToC(const arma::mat& indsR) {
    arma::umat indsC(indsR.n_rows, indsR.n_cols);
    std::transform(indsR.begin(), indsR.end(), indsC.begin(), [](double x) { return static_cast<std::size_t>(x) - 1; });
    return indsC;
}

arma::vec computeLevelSizes(const std::size_t levelCount, const arma::uvec& inds) {
    arma::vec levelSizes(levelCount, arma::fill::zeros);
    for (auto x : inds)
        levelSizes[x] += 1.0;
    return levelSizes;
}

std::unique_ptr<const FixedEffects> FixedEffects::create(const arma::uvec& levelCounts, const arma::mat& indsR) {
    auto effects = std::make_unique<FixedEffects>();
    effects->size = levelCounts.size();
    effects->indicators = createIndicators(levelCounts, indsR);

    effects->simpleEffects = std::vector<SimpleFixedEffect>();
    for (auto i = 0u; i < effects->size; i ++)
        effects->simpleEffects.emplace_back(effects->indicators[i], std::move(singleton));

    auto indsC = convertNumberingRToC(indsR);
    effects->componentTables = computeComponents(levelCounts, indsC);

    return effects;
}

std::vector<CrossComponentError> FixedEffects::checkComponents(const arma::mat& indsR) const {
    auto indsC = convertNumberingRToC(indsR);
    if (!componentTables)
        return {};

    return ::checkComponents(*componentTables, indsC);
}

struct Payload {
    const FixedEffects& fixedEffects;
    arma::subview_col<double> data;
    arma::vec backup;
    std::vector<arma::vec> deltas;
    const double epsilon;
    const std::size_t maxIterations;

    Payload(const FixedEffects& fixedEffects, arma::subview_col<double> data, const double epsilon, const std::size_t maxIterations): fixedEffects(fixedEffects), data(data), backup(data.n_cols), deltas(), epsilon(epsilon), maxIterations(maxIterations) {
        std::transform(fixedEffects.simpleEffects.begin(),
                       fixedEffects.simpleEffects.end(),
                       std::back_inserter(deltas),
                       [](auto& x) { return arma::zeros(x.indicator.levelCount); });
    }
};

void demean(void* ptr) {
    auto& casted = *(static_cast<Payload*>(ptr));

    auto i = 0;
    do {
        casted.backup = casted.data;
        for (auto i = 0u; i < casted.fixedEffects.size; i ++)
            casted.deltas[i] += casted.fixedEffects.simpleEffects[i].demean(casted.data);
        if (arma::accu(arma::abs(casted.data - casted.backup)) < casted.epsilon)
            return;
        i ++;
    } while (i < casted.maxIterations);

    mainQueue->commit(ptr);
}

std::vector<arma::mat> FixedEffects::demean(arma::mat& data) const {
    std::vector<std::unique_ptr<Payload>> payloads;
    for (auto i = 0u; i < data.n_cols; i ++) {
        auto payload = std::make_unique<Payload>(*this, data.col(i), 1e-6, 5);
        payloads.emplace_back(std::move(payload));
    }

    std::queue<void*> queue;
    for (const auto& payload : payloads)
        queue.push(static_cast<void*>(payload.get()));
    mainQueue->commit(::demean, std::move(queue));
    mainQueue->crush();

    std::vector<arma::mat> deltas;
    for (auto i = 0u; i < size; i ++) {
        const auto& effect = simpleEffects[i];
        arma::mat delta(effect.indicator.levelCount, data.n_cols);
        for (auto j = 0u; j < data.n_cols; j ++)
            delta.col(j) = payloads[j]->deltas[i];
        deltas.emplace_back(std::move(delta));
    }
    return deltas;
}
