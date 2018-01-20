#include "FastFESolver.h"

FastFESolver::FastFESolver(const arma::mat& X, const arma::colvec& Y,
                           const std::vector<FixedEffect>& fixedEffects):
    PlainModel(X, Y, fixedEffects) {}

void FastFESolver::compute() {
    demean();
    arma::mat initData = arma::join_rows(Y, X);
    arma::mat projData = arma::join_rows(Y_, X_);
    auto result = estimate(initData, projData, deltas);
    this->result = PlainModel::Result();
    this->result.params = result.projModel.beta;
    this->result.effects = result.fixedEffects;
    this->result.residuals = result.residuals;
    this->result.fittedValues = result.fittedValues;
    this->result.intercept = result.intercept;
}

void FastFESolver::demean() {
    arma::mat data = arma::join_rows(Y, X);
    arma::mat copy;

    for (const auto& category : fixedEffects)
        deltas.push_back(arma::zeros(category.groupCount, data.n_cols));

    do {
        copy = data;
        for (int i = 0; i < fixedEffects.size(); i ++)
            fixedEffects[i].demean(data, deltas[i]);
    } while (arma::accu(abs(data - copy)) > 1e-5);

    Y_ = data.col(0);
    if (data.n_cols > 1)
      X_ = data.cols(1, data.n_cols - 1);
    else
      X_ = arma::mat(data.n_rows, 0);
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

Result estimate(const arma::mat& initData, const arma::mat& projData, const std::vector<arma::mat>& deltas) {
    auto projModel = LinearModel::solve(projData);

    double intercept = 0;
    auto _1 = estimateFixedEffects(deltas, projModel.beta);
    auto fixedEffects = _1.first;
    intercept += _1.second;
    auto _2 = estimateResiduals(projData, projModel.beta);
    auto residuals = _2.first;
    intercept += _2.second;
    auto fittedValues = estimateFittedValues(initData, residuals);

    return Result { projModel, fixedEffects, residuals, fittedValues, intercept };
}
