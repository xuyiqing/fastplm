#include "FixedEffect.h"
#include "GradientDescent.h"

// [[Rcpp::export()]]
List fastplm(arma::mat rawData, arma::mat rawFixedEffects) {
  arma::mat X = rawData.cols(1, rawData.n_cols - 1);
  arma::mat Y = rawData.col(0);
  
  std::vector<FixedEffect> fixedEffects;
  for (int i = 0; i < rawFixedEffects.n_cols; i ++)
    fixedEffects.push_back(FixedEffect::fromColumn(rawFixedEffects.col(i)));

  GradientDescent solver(X, Y, fixedEffects);
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
