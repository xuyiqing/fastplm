#include "FixedEffectModel.h"

auto estimateFixedEffects(const std::vector<FixedEffects::Deltas>& deltas, const arma::vec& coefficients) {
    double totalIntercept = 0;

    std::vector<arma::vec> sfeCoefs;
    auto simpleEffectCount = deltas[0].first.size();
    for (auto i = 0u; i < simpleEffectCount; i ++) {
        arma::vec deltaY = deltas[0].first[i];

        for (auto j = 0u; j < coefficients.n_elem; j ++)
            deltaY -= coefficients[j] * deltas[j + 1].first[i];

        double intercept = arma::mean(deltaY);
        deltaY -= intercept;
        totalIntercept += intercept;
        sfeCoefs.push_back(std::move(deltaY));
    }

    std::vector<arma::mat> cfeCoefs;
    auto complexEffectCount = deltas[0].second.size();
    for (auto i = 0u; i < complexEffectCount; i ++) {
        arma::mat deltaY = deltas[0].second[i];

        for (auto j = 0; j < coefficients.n_elem; j ++)
            deltaY -= coefficients[j] * deltas[j + 1].second[i];

        cfeCoefs.push_back(deltaY.t());
    }

    return std::make_tuple(sfeCoefs, cfeCoefs, totalIntercept);
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

const FixedEffectModel FixedEffectModel::solve(const arma::mat& initData, const FixedEffects& fixedEffects) {
    FixedEffectModel model;

    arma::mat projData = initData;
    auto deltas = fixedEffects.demean(projData);
    model.demeaned = LinearModel::solve(projData);

    model.intercept = 0;
    auto _1 = estimateFixedEffects(deltas, model.demeaned.beta);
    model.sfeCoefs = std::get<0>(_1);
    model.cfeCoefs = std::get<1>(_1);
    model.intercept += std::get<2>(_1);
    auto _2 = estimateResiduals(projData, model.demeaned.beta);
    model.residuals = _2.first;
    model.intercept += _2.second;
    model.fittedValues = estimateFittedValues(initData, model.residuals);

    return model;
}
