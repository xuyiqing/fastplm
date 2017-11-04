#include "CrashQueue.h"
#include "FastFESolver.h"
#include "SlowFESolver.h"
#include "GPSolver.h"

// [[Rcpp::export()]]
List internalSolveFE(arma::mat rawData, arma::mat rawFixedEffects, unsigned coreNum = 1, bool estimateFE = false) {
  mainQueue = new CrashQueue(coreNum);
  ScopeGuard _([]{ delete mainQueue; mainQueue = nullptr; });
  
  arma::mat X = rawData.cols(1, rawData.n_cols - 1);
  arma::mat Y = rawData.col(0);
  
  std::vector<FixedEffect> fixedEffects;
  for (int i = 0; i < rawFixedEffects.n_cols; i ++)
    fixedEffects.push_back(FixedEffect::fromColumn(rawFixedEffects.col(i)));
    
  FastFESolver solver(X, Y, fixedEffects, estimateFE);
  solver.compute();
  
  List result;
  result["coefficients"] = solver.result.params;
  result["fittedValues"] = solver.result.fittedValues;
  result["residuals"] = solver.result.residuals;
  
  if (estimateFE) {
    result["intercept"] = solver.result.intercept;
    
    List FEcoefs;
    for (int i = 0; i < solver.result.effects.size(); i ++) {
      auto effectId = "effect" + std::to_string(i + 1);
      FEcoefs[effectId] = solver.result.effects[i];
    }
    result["FEcoefs"] = FEcoefs;
  }
  
  return result;
}

// [[Rcpp::export()]]
List solveGP(arma::vec Y, arma::mat X, arma::mat tois, arma::mat iots, bool isBalanced) {
    auto timeCount = tois.n_rows, indivCount = iots.n_rows;
    arma::mat Y_(Y.memptr(), timeCount, indivCount, false);
    arma::cube X_(X.memptr(), timeCount, indivCount, X.n_cols, false);
    GPSolver solver(X_, Y_, tois, iots, isBalanced);
    
    //Rcout << X_ << "\n";
    Rcout << Y_ << "\n";
    
    List result;
    result["coefficients"] = solver.compute();
    return result;
}


