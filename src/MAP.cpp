#include "MAP.h"
#include "GradientDescent.h"

MAP::MAP(const arma::mat& X, const arma::colvec& Y,
         const std::vector<FixedEffect>& fixedEffects,
         DemeanTransform transform):
    PlainModel(X, Y, fixedEffects), transform(transform) {}

void MAP::compute() {
    demean();
    estimateParams();
    estimateFixedEffects();
}

void MAP::demean() {
    arma::mat data = arma::join_rows(Y, X);
    arma::mat copy;
    
    do {
        copy = data;
        switch (transform) {
            case kHalperin:
                halperin(data);
                break;
            case kSymmetricHalperin:
                symmetricHalperin(data);
                break;
            case kCimmino:
                throw std::logic_error("Cimmino not ready yet");
                break;
            default:
                throw std::logic_error("Unknown transform");
        }
    } while (arma::accu(abs(data - copy)) > 1e-5);
    
    Y_ = data.col(0);
    X_ = data.cols(1, data.n_cols - 1);
}

void MAP::halperin(arma::mat& data) {
    for (const auto& category: fixedEffects)
        demean(data, category);
}

void MAP::symmetricHalperin(arma::mat& data) {
    int i = 0;
    for (; i < fixedEffects.size(); i ++)
        demean(data, fixedEffects[i]);
    for (; i >= 0; i --)
        demean(data, fixedEffects[i]);
}

void MAP::demean(arma::mat& data, const FixedEffect& category) {
    for (int i = 0; i < category.groupCount; i ++) {
        auto indicator = category.indicator.col(i);
        int count = category.valuesOccurences[i];
        if (count > 0) {
            auto means = (indicator.t() * data) / static_cast<double>(count);
            data -= (indicator * means);
        }
    }
}

void MAP::estimateParams() {
    // Remove columns of X_ that has no variation.
    std::vector<int> badCols;
    arma::mat goodCols(X_.n_rows, 0);
    for (int i = 0; i < X_.n_cols; i ++) {
        arma::colvec var = arma::unique(X_.col(i));
        if (var.n_rows == 1)
            badCols.push_back(i);
        else
            goodCols.insert_cols(goodCols.n_cols, X_.col(i));
    }
    
    const arma::colvec goodParams = solve(goodCols, Y);
    
    if (badCols.size() == 0) {
        result.params = goodParams;
        return;
    }
    
    result.params = arma::colvec(X.n_cols, 1);
    for (int i = 0, iBad = 0, iGood = 0; i < paramCount; i ++) {
        if (i == badCols[iBad]) {
            result.params(i) = arma::datum::nan;
            iBad ++;
        } else {
            result.params(i) = goodParams(iGood ++);
        }
    }
}

void MAP::estimateFixedEffects() {
    if (fixedEffects.size() == 0) {
        result.effects = std::vector<arma::colvec>();
        return;
    }
    
    auto residuals = Y - X * result.params;
    GradientDescent solver(arma::zeros(observationCount, paramCount), residuals, fixedEffects, false);
    solver.compute();
    result.effects = solver.result.effects;
    result.intercept = solver.result.intercept;
}
