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
    this->result.params = result.beta;
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

const double DOUBLE_TOLERANCE = 1e-5;

inline bool isZero(const arma::colvec& vec, const double eps = DOUBLE_TOLERANCE) {
    return arma::all(arma::abs(vec) < eps);
}

arma::mat orthogonalize(arma::mat mat) {
    for (auto i = 0u; i < mat.n_cols; i ++) {
        arma::colvec v = mat.col(i);

        for (auto j = 0u; j < i; j ++)
            v -= arma::dot(mat.col(j), mat.col(i)) * mat.col(j);

        mat.col(i) = v / arma::norm(v);
    }
    return mat;
}

auto checkLinearDependency(const arma::mat& mat) {
    auto ortho = orthogonalize(mat);

    arma::colvec dependents = arma::zeros(mat.n_cols);
    arma::colvec independents = arma::zeros(mat.n_cols);

    for (auto i = 0u; i < mat.n_cols; i ++)
        if (isZero(ortho.col(i)))
            dependents[i] = 1;
        else
            independents[i] = 1;

    return std::make_pair(arma::find(dependents), arma::find(independents));
}

inline const arma::colvec getY(const arma::mat& data) {
    return data.col(0);
}

inline const arma::mat getX(const arma::mat& data) {
    return data.cols(1, data.n_cols - 1);
}

arma::colvec estimateBeta(const arma::mat& data, const arma::uvec& independents) {
    arma::colvec leanBeta = solve(getX(data).cols(independents), getY(data));
    arma::colvec fatBeta = arma::zeros(data.n_cols - 1);
    fatBeta.rows(independents) = leanBeta;
    return fatBeta;
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
    auto tmp0 = checkLinearDependency(getX(projData));
    arma::uvec dependents = tmp0.first;
    arma::uvec independents = tmp0.second;

    double intercept = 0;
    auto beta = estimateBeta(projData, independents);
    auto tmp1 = estimateFixedEffects(deltas, beta);
    auto fixedEffects = tmp1.first;
    intercept += tmp1.second;
    auto tmp2 = estimateResiduals(projData, beta);
    auto residuals = tmp2.first;
    intercept += tmp2.second;
    auto fittedValues = estimateFittedValues(initData, residuals);

    return Result { beta, dependents, independents, fixedEffects, residuals, fittedValues, intercept };
}
