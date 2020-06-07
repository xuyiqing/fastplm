#include <queue>

#include "CrushQueue.h"
#include "FixedEffects.h"
#include <cmath>
using std::isnan;

std::unique_ptr<const FixedEffects> FixedEffects::create(const arma::uvec& levelCounts, const arma::mat& indsR, const arma::uvec& simpleEffects, const arma::uvec& complexEffects, const arma::uvec& complexInfluences, const std::vector<arma::mat>& weights) {
    auto effects = std::make_unique<FixedEffects>();

    effects->indicators = createIndicators(levelCounts, indsR);
    effects->simpleEffects = std::vector<SimpleFixedEffect>();
    for (auto i = 0u; i < simpleEffects.n_elem; i ++) {
        const auto& eff = effects->indicators[simpleEffects[i]];
        effects->simpleEffects.emplace_back(eff, SimpleInfluence());
    }

    const auto totalRows = indsR.n_rows;

    for (auto i = 0u; i < complexEffects.n_elem; i ++) {
        const auto& eff = effects->indicators[complexEffects[i]];
        const auto& inf = effects->indicators[complexInfluences[i]];
        arma::mat weight = weights[i];
        const auto infDim = weight.n_rows;

        std::vector<arma::mat> normalizers;
        for (auto effId = 0u; effId < eff.levelCount; effId ++) {
            const auto levelSize = eff.levelSizes[effId];
            arma::mat projection(infDim, levelSize);
            auto col = 0;

            for (auto row = 0u; row < totalRows; row ++) {
                if (eff.indicator[row] != effId)
                    continue;

                auto infId = inf.indicator[row];
                projection.col(col ++) = weight.col(infId);
            }

            const arma::mat normalizer = arma::inv(projection * projection.t());
            normalizers.push_back(std::move(normalizer));
        }

        ComplexInfluence CI(inf, std::move(weight), std::move(normalizers));
        effects->complexEffects.emplace_back(eff, std::move(CI));
    }

    effects->componentTables = computeComponents(levelCounts, effects->indicators);

    return effects;
}

std::vector<CrossComponentError> FixedEffects::checkComponents(const arma::mat& indsR) const {
    arma::umat indsC(indsR.n_rows, indsR.n_cols);
    std::transform(indsR.begin(), indsR.end(), indsC.begin(), [](double x) { return static_cast<std::size_t>(x) - 1; });

    if (!componentTables)
        return {};

    return ::checkComponents(*componentTables, indsC);
}

struct Payload {
    const FixedEffects& fixedEffects;
    arma::subview_col<double> data;
    arma::vec backup;
    std::vector<arma::vec> simpleDeltas;
    std::vector<arma::mat> complexDeltas;
    const double epsilon;
    const std::size_t maxIterations;

    Payload(const FixedEffects& fixedEffects, arma::subview_col<double> data, const double epsilon, const std::size_t maxIterations): fixedEffects(fixedEffects), data(data), backup(data.n_cols), simpleDeltas(), epsilon(epsilon), maxIterations(maxIterations) {
        std::transform(fixedEffects.simpleEffects.begin(),
                       fixedEffects.simpleEffects.end(),
                       std::back_inserter(simpleDeltas),
                       [](auto& x) { return arma::zeros(x.indicator.levelCount); });

        std::transform(fixedEffects.complexEffects.begin(),
                       fixedEffects.complexEffects.end(),
                       std::back_inserter(complexDeltas),
                       [](auto& x) { return arma::zeros(x.influence.dimension, x.indicator.levelCount); });
    }
};

void demean(void* ptr) {
    auto& casted = *(static_cast<Payload*>(ptr));
    const auto& effects = casted.fixedEffects;

    auto i = 0;
    do {
        casted.backup = casted.data;

        for (auto i = 0u; i < effects.simpleEffects.size(); i ++)
            casted.simpleDeltas[i] += effects.simpleEffects[i].demean(casted.data);
        for (auto i = 0u; i < effects.complexEffects.size(); i ++)
            casted.complexDeltas[i] += effects.complexEffects[i].demean(casted.data);

        double epsilon = arma::accu(arma::abs(casted.data - casted.backup));
        if (isnan(epsilon))
            //exit(-1);
            Rcpp::stop("Demean Fails.\n");
        if (epsilon < casted.epsilon)
            return;
        i ++;
    } while (i < casted.maxIterations);

    mainQueue->commit(ptr);
}

std::vector<FixedEffects::Deltas> FixedEffects::demean(arma::mat& data) const {
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

    std::vector<Deltas> deltas;
    for (auto i = 0u; i < data.n_cols; i ++)
        deltas.emplace_back(std::make_pair(std::move(payloads[i]->simpleDeltas),
                                           std::move(payloads[i]->complexDeltas)));

    return deltas;
}
