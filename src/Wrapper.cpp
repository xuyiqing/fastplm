#include "CrushQueue.h"
#include "FixedEffectModel.h"
#include "GPSolver.h"

using Rcpp::List;

// [[Rcpp::export()]]
List internalSolveFE(arma::mat rawData, arma::mat rawFixedEffects, unsigned coreNum = 1) {
  mainQueue = new CrushQueue(coreNum);
  ScopeGuard _([]{ delete mainQueue; mainQueue = nullptr; });
  
  arma::mat X;
  if (rawData.n_cols > 1)
    X = rawData.cols(1, rawData.n_cols - 1);
  else
    X = arma::mat(rawData.n_rows, 0);
  arma::mat Y = rawData.col(0);
  
  std::vector<FixedEffect> fixedEffects;
  for (int i = 0; i < rawFixedEffects.n_cols; i ++)
    fixedEffects.push_back(FixedEffect::fromColumn(rawFixedEffects.col(i)));

  auto result = FixedEffectModel::solve(rawData, fixedEffects);
  return static_cast<List>(result);
  return result;
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


