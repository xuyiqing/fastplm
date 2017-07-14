#include "FastFESolver.h"
#include "SlowFESolver.h"
#include "GPSolver.h"

// [[Rcpp::export()]]
List solveFE(arma::mat rawData, arma::mat rawFixedEffects) {
  arma::mat X = rawData.cols(1, rawData.n_cols - 1);
  arma::mat Y = rawData.col(0);
  
  std::vector<FixedEffect> fixedEffects;
  for (int i = 0; i < rawFixedEffects.n_cols; i ++)
    fixedEffects.push_back(FixedEffect::fromColumn(rawFixedEffects.col(i)));
    
  FastFESolver solver(X, Y, fixedEffects);
  solver.compute();
  
  List result;
  result["coefficients"] = solver.result.params;
  result["intercept"] = solver.result.intercept;
  for (int i = 0; i < solver.result.effects.size(); i ++) {
    auto effectId = "effect" + std::to_string(i + 1);
    result[effectId] = solver.result.effects[i];
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


