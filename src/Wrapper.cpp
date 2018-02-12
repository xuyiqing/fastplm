#include "CrushQueue.h"
#include "FixedEffects.h"
#include "FixedEffectModel.h"
#include "GPSolver.h"

using Rcpp::List;
using Rcpp::XPtr;

// [[Rcpp::export()]]
SEXP CreateFixedEffects(arma::uvec groupSizes, arma::mat indicators) {
    auto ptr = FixedEffects::create(groupSizes, indicators);
    return XPtr<const FixedEffects>(ptr.release());
}

// [[Rcpp::export()]]
bool ContainMultipleComponents(SEXP wrappedFixedEffects) {
    auto fixedEffects = XPtr<const FixedEffects>(wrappedFixedEffects);
    return static_cast<bool>(fixedEffects->componentTables);
}

// [[Rcpp::export()]]
List CheckComponents(SEXP wrappedFixedEffects, arma::mat indicators) {
    auto fixedEffects = XPtr<const FixedEffects>(wrappedFixedEffects);
    auto cppErrors = fixedEffects->checkComponents(indicators);

    List rErrors;
    for (const auto& error : cppErrors)
        rErrors.push_back(static_cast<List>(error));
    return rErrors;
}


// [[Rcpp::export()]]
List SolveFixedEffects(arma::mat data, SEXP wrappedFixedEffects, std::size_t coreNum = 1) {
    mainQueue = new CrushQueue(coreNum);
    ScopeGuard _([]{ delete mainQueue; mainQueue = nullptr; });

    auto fixedEffects = XPtr<const FixedEffects>(wrappedFixedEffects);
    return FixedEffectModel::solve(data, *fixedEffects);
}

// [[Rcpp::export()]]
List solveGP(arma::vec Y, arma::mat X, arma::mat tois, arma::mat iots, bool isBalanced) {
    auto timeCount = tois.n_rows, indivCount = iots.n_rows;
    arma::mat Y_(Y.memptr(), timeCount, indivCount, false);
    arma::cube X_(X.memptr(), timeCount, indivCount, X.n_cols, false);
    GPSolver solver(X_, Y_, tois, iots, isBalanced);
    
    List result;
    result["coefficients"] = solver.compute();
    return result;
}
