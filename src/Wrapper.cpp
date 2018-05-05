#include "CrushQueue.h"
#include "FixedEffects.h"
#include "FixedEffectModel.h"

using Rcpp::List;
using Rcpp::XPtr;

// [[Rcpp::export()]]
SEXP CreateFixedEffects(arma::uvec levelCounts, arma::mat indsR,
                        arma::uvec simpleEffects,
                        arma::uvec complexEffects,
                        arma::uvec complexInfluences,
                        List wrappedWeights) {
    simpleEffects -= 1u;
    complexEffects -= 1u;
    complexInfluences -= 1u;

    std::vector<arma::mat> weights;
    for (const SEXP& elem : wrappedWeights) {
        arma::mat m = Rcpp::as<arma::mat>(elem);
        weights.push_back(std::move(m));
    }

    auto ptr = FixedEffects::create(levelCounts, indsR,
        simpleEffects, complexEffects, complexInfluences, weights);

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
