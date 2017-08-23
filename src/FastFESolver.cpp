#include "FastFESolver.h"
#include "SlowFESolver.h"

FastFESolver::FastFESolver(const arma::mat& X, const arma::colvec& Y,
                           const std::vector<FixedEffect>& fixedEffects,
                           bool doesComputeFixedEffects):
    PlainModel(X, Y, fixedEffects), doesComputeFixedEffects(doesComputeFixedEffects){}

void FastFESolver::compute() {
    demean();
    estimateParams();
    if (doesComputeFixedEffects)
        estimateFixedEffects();
}

void FastFESolver::demean() {
    arma::mat data = arma::join_rows(Y, X);
    arma::mat copy;
    
    do {
        copy = data;
        for (const auto& category: fixedEffects)
            category.demean(data);
    } while (arma::accu(abs(data - copy)) > 1e-5);
    
    Y_ = data.col(0);
    X_ = data.cols(1, data.n_cols - 1);
}

void FastFESolver::estimateParams() {
    // Remove columns of X_ that has no variation.
    std::vector<int> badCols;
    
    for (int i = 0; i < X_.n_cols; i ++) {
        bool noVar = true;
        double bookmark = X_(0, i);
        for (int j = 0; j < X_.n_rows; j ++)
            if (bookmark != X_(j, i)) {
                noVar = false;
                break;
            }
        
        if (noVar)
            badCols.push_back(i);
    }
    
    for (int i = badCols.size() - 1; i >= 0; i --)
        X_.shed_col(i);
    
    const arma::colvec params = solve(X_, Y_);
    
    if (badCols.size() == 0) {
        result.params = params;
    }
    else {
        result.params = arma::colvec(X.n_cols, 1);
        for (int i = 0, iBad = 0, iGood = 0; i < paramCount; i ++) {
            if (i == badCols[iBad]) {
                result.params(i) = arma::datum::nan;
                iBad ++;
            } else {
                result.params(i) = params(iGood ++);
            }
        }
    }
    
    result.residuals = Y_ - X_ * result.params;
    result.fittedValues = Y - result.residuals;
}

void FastFESolver::estimateFixedEffects() {
    if (fixedEffects.size() == 0) {
        result.effects = std::vector<arma::colvec>();
        return;
    }
    
    auto residuals = Y - X * result.params;
    SlowFESolver solver(arma::zeros(observationCount, paramCount), residuals, fixedEffects, false);
    solver.compute();
    result.effects = solver.result.effects;
    result.intercept = solver.result.intercept;
}
