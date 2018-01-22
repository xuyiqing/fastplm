#include "FixedEffectModel.h"

auto demean(const arma::mat& initData, const std::vector<FixedEffect>& fixedEffects) {
    arma::mat data = initData;
    arma::mat copy;
    std::vector<arma::mat> deltas;

    for (const auto& category : fixedEffects)
        deltas.push_back(arma::zeros(category.groupCount, data.n_cols));

    do {
        copy = data;
        for (int i = 0; i < fixedEffects.size(); i ++)
            fixedEffects[i].demean(data, deltas[i]);
    } while (arma::accu(abs(data - copy)) > 1e-5);

    return std::make_pair(data, deltas);
}

auto estimateFixedEffects(const std::vector<arma::mat>& deltas, const arma::colvec& beta) {
    std::vector<arma::colvec> fixedEffects;
    fixedEffects.reserve(deltas.size());

    double totalIntercept = 0.0;

    std::transform(deltas.cbegin(), deltas.cend(), std::back_inserter(fixedEffects),
                   [&](const auto& delta) -> arma::colvec {
                       arma::colvec fixedEffect = getY(delta) - getX(delta) * beta;
                       double intercept = arma::mean(fixedEffect);
                       fixedEffect -= intercept;
                       totalIntercept += intercept;
                       return fixedEffect;
                   });

    return std::make_pair(fixedEffects, totalIntercept);
}

inline auto estimateResiduals(const arma::mat& data, const arma::colvec& beta) {
    arma::colvec residuals = getY(data) - getX(data) * beta;
    double intercept = arma::mean(residuals);
    residuals -= intercept;
    return std::make_pair(residuals, intercept);
}

inline arma::colvec estimateFittedValues(const arma::mat& initData, const arma::colvec& residuals) {
    return getY(initData) - residuals;
}

const FixedEffectModel FixedEffectModel::solve(const arma::mat& initData, const std::vector<FixedEffect>& fixedEffects) {
    FixedEffectModel model;

    auto _0 = demean(initData, fixedEffects);
    auto projData = _0.first;
    auto deltas = _0.second;
    model.demeaned = LinearModel::solve(projData);

    model.intercept = 0;
    auto _1 = estimateFixedEffects(deltas, model.demeaned.beta);
    model.feCoefs = _1.first;
    model.intercept += _1.second;
    auto _2 = estimateResiduals(projData, model.demeaned.beta);
    model.residuals = _2.first;
    model.intercept += _2.second;
    model.fittedValues = estimateFittedValues(initData, model.residuals);

    return model;
}
