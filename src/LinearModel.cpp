#include "LinearModel.h"

arma::mat orthogonalize(arma::mat mat) {
    for (auto i = 0u; i < mat.n_cols; i ++) {
        arma::vec v = mat.col(i);

        for (auto j = 0u; j < i; j ++)
            v -= arma::dot(mat.col(j), mat.col(i)) * mat.col(j);

        mat.col(i) = v / arma::norm(v);
    }
    return mat;
}

auto checkLinearDependency(const arma::mat& mat) {
    auto ortho = orthogonalize(mat);

    arma::vec dependents = arma::zeros(mat.n_cols);
    arma::vec independents = arma::zeros(mat.n_cols);

    for (auto i = 0u; i < mat.n_cols; i ++)
        if (isZero(ortho.col(i)))
            dependents[i] = 1;
        else
            independents[i] = 1;

    return std::make_pair(arma::find(dependents), arma::find(independents));
}

const LinearModel LinearModel::solve(const arma::mat& X, const arma::vec& Y) {
    auto _1 = checkLinearDependency(X);
    arma::uvec dependents = _1.first;
    arma::uvec independents = _1.second;
    bool isLinearDependent = dependents.n_rows > 0;

    arma::vec leanBeta = arma::solve(X.cols(independents), Y);
    arma::vec fatBeta = arma::zeros(X.n_cols);
    fatBeta.rows(independents) = leanBeta;

    LinearModel model;
    model.X = X;
    model.Y = Y;
    model.isLinearDependent = isLinearDependent;
    model.dependents = dependents;
    model.independents = independents;
    model.beta = fatBeta;
    return model;
}

const LinearModel LinearModel::solve(const arma::mat& data) {
    if (data.n_cols > 1)
        return solve(getX(data), getY(data));

    LinearModel model;
    model.X = getX(data);
    model.Y = getY(data);
    model.isLinearDependent = 0;

    model.dependents = model.dependents = arma::zeros<arma::uvec>(0);
    model.beta = arma::zeros<arma::vec>(0);
    return model;
}
