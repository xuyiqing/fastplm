// [[Rcpp::plugins(cpp11)]]

#include <vector>
#include <algorithm>
#include <unordered_map>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

#include "Demean.h"

/*
*  Terminology
*
*  Data Matrix
*  A data matrix consist of many *observations*, each as a row vector;
*  A data matrix consist of a column of results (Y), and many columns of variables (X_1, ..., X_N)
*  observation[0] is *the* observation, the result of some observation;
*  observation[i] records the *i-th variable* of that observation.
*
*  Fixed Effects
*  A matrix of fixed effects consist of several *categories*, each as a column vector;
*  category[i] records the *fixed effect value* of this category for the i-th observation.
*/

std::tuple<arma::colvec, arma::colvec, arma::colvec>
estimateCovariates(const arma::colvec &Y, const arma::mat &X, const int groupCount) {
  // Remove columns of X that has no variation.
  std::vector<int> badCols;
  arma::mat goodX(X.n_rows, 0);
  for (int i = 0; i < X.n_cols; i ++) {
    arma::colvec var = arma::unique(X.col(i));
    if (var.n_rows == 1)
      badCols.push_back(i);
    else
      goodX.insert_cols(goodX.n_cols, X.col(i));
  }
  
  const int covariateCount = goodX.n_cols;
  const arma::colvec goodCoefficients = solve(goodX, Y);
  const arma::colvec residuals = Y - goodX * goodCoefficients;
  
  const int degreeOfFreedom = X.n_rows - groupCount - covariateCount;
  const double sigma2 = arma::as_scalar(residuals.t() * residuals) / static_cast<double>(degreeOfFreedom);
  const arma::colvec goodStdErrors = arma::sqrt(sigma2 * arma::diagvec(arma::inv(X.t() * X)));
  
  if (badCols.size() == 0)
    return std::make_tuple(residuals, goodCoefficients, goodStdErrors);
  
  arma::colvec coefficients(X.n_cols, 1);
  arma::colvec stdErrors(X.n_cols, 1);
  for (int i = 0, iGood = 0, iBad = 0;
       i < X.n_cols;
       i ++) {
    if (i == badCols[iBad]) {
      coefficients(i) = arma::datum::nan;
      stdErrors(i) = arma::datum::nan;
      iBad ++;
    } else {
      coefficients[i] = goodCoefficients(iGood);
      stdErrors(i) = goodStdErrors(iGood);
      iGood ++;
    }
  }
  
  return std::make_tuple(residuals, coefficients, stdErrors);
}

// Not a strict Laplacian.
std::pair<arma::colvec, arma::mat> buildLaplacian(const arma::colvec &residuals,
                                                  const arma::colvec &cat1,
                                                  const arma::colvec &cat2) {
  std::unordered_map<double, int> catIndices1, catIndices2;
  int ix = 0;
  for (auto group : cat1) {
    if (catIndices1.count(group) == 0)
      catIndices1[group] = ix ++;
  }
  for (auto group : cat2) {
    if (catIndices2.count(group) == 0)
      catIndices2[group] = ix ++;
  }
  int groupCount = ix;
  arma::colvec LY = arma::zeros(groupCount);
  arma::mat LX = arma::zeros(groupCount, groupCount);
  
  for (int i = 0; i < cat1.n_rows; i ++) {
    int group1Index = catIndices1[cat1(i)];
    int group2Index = catIndices2[cat2(i)];
    double y = residuals(i);
    LY(group1Index) += y;
    LY(group2Index) -= y;
    LX(group1Index, group1Index) ++;
    LX(group1Index, group2Index) --;
    LX(group2Index, group2Index) ++;
    LX(group2Index, group1Index) --;
  }
  
  return std::make_pair(LY, LX);
}

// [[Rcpp::export()]]
arma::colvec playAround(const arma::colvec R, const arma::colvec C1, const arma::colvec C2) {
  auto system = buildLaplacian(R, C1, C2);
  return solve(system.second, system.first);
}

// [[Rcpp::export()]]
List wrapper(const arma::colvec R, const arma::colvec C1, const arma::colvec C2) {
  auto system = buildLaplacian(R, C1, C2);
  List output;
  output["X"] = system.second;
  output["Y"] = system.first;
  return output;
}
arma::colvec computeFixedEffects(const arma::mat &data, arma::colvec coefficients, const arma::mat &fixedEffects) {
  auto Y = data.col(0);
  auto X = data.cols(1, data.n_cols - 1);
  arma::colvec residualWithFEs = Y - X * coefficients;
  auto Laplacian = buildLaplacian(residualWithFEs, fixedEffects.col(0), fixedEffects.col(1));
  arma::colvec alphas = solve(Laplacian.second, Laplacian.first);
  return alphas;
}

// [[Rcpp::export()]]
List fastplm(arma::mat data,
             arma::mat FE,
             int FEcoefs = 0
){
  
  // parse data
  int n = data.n_rows;
  int k = data.n_cols;
  int m = FE.n_cols;
  int p = k-1; // No. of covariates
  arma::mat data_bak = data;
  
  
  // count total number of groups (loss of degrees of freedom)
  arma::mat FEvalues;
  for(int ii=0; ii<m; ii++){ // cluster
    arma::colvec fe = arma::unique(FE.col(ii));
    int g=fe.n_rows; // No. of group values
    arma::colvec gp(g);
    gp.fill(ii);
    arma::mat FEvalue = join_rows(gp, fe);
    FEvalues = join_cols(FEvalues, FEvalue);
  }
  int gtot = FEvalues.n_rows;
  
  data = demeanMAP(data, FE);
  
  // declare variables
  arma::colvec resid; // n*1
  arma::colvec coeff; // coefficient (full)
  arma::colvec se; // SE (full)
  
  if (p > 0) {
    auto estimationResult = estimateCovariates(data.col(0), data.cols(1, p), gtot);
    resid = std::get<0>(estimationResult);
    coeff = std::get<1>(estimationResult);
    se = std::get<2>(estimationResult);
  } else {
    resid = data.col(0);
  }
  
  // auto alphas = computeFixedEffects(data_bak.col(0), data_bak.cols(1, p), coeff, FE);
  arma::colvec y;  // n*1
  arma::mat X; // n*p
  arma::colvec e; // n*1 (with fixed effects)
  arma::colvec coef; // coefficient
  double mu; // grand mean
  arma::colvec LHS; // group means
  arma::mat W; // big weighting matrix to calculate fixed effects
  arma::colvec alphas; // fixed effect coefficients
  
  
  // Calculate fixed effects coefficients
  if (FEcoefs == 1) {
    data = data_bak;
    y = data.col(0);  // n*1
    
    // grand mean
    mu = arma::mean(y);
    if (p > 0) {
      X = data.cols(1, p); // n*p
      coef = coeff;
      for (int i=0; i<p; i++) {
        if (coef(i) == arma::datum::nan) {
          coef(i) = 0;
        }
      }
      mu = mu - arma::as_scalar(arma::mean(X, 0) * coef);
    }
    
    // residuals (with fixed effects)
    e = y - mu;
    if (p > 0) {
      e = e - X * coef;
    }
    
    arma::colvec LHS(gtot, arma::fill::zeros);
    arma::mat W(gtot, gtot, arma::fill::zeros);
    arma::colvec FEval1 = arma::unique(FE.col(0));
    arma::colvec FEval2 = arma::unique(FE.col(1));
    int f1 = FEval1.n_rows;
    int f2 = FEval2.n_rows;
    
    for (int i=0; i<n; i++) {
      int cont = 1; // continue
      int j = 0;
      while ((cont==1) & (j < f1)) {
        if (FE(i,0) == FEval1[j]) {
          int k = 0;
          while ((cont==1) & (k < f2)) {
            if (FE(i,1) == FEval2[k]) {
              LHS(j) = LHS(j) + e(i);
              LHS(f1 + k) = LHS(f1 + k) + e(i);
              W(j, j) = W(j, j) + 1;
              W(j, f1 + k) = W(j, f1 + k) + 1;
              W(f1 + k, j) = W(f1 + k, j) + 1;
              W(f1 + k, f1 +k) = W(f1 + k, f1 +k) + 1;
              cont = 0;
            }
            k = k + 1;
          }
        }
        j = j + 1;
      }
    }
    alphas = arma::solve(W,LHS);
    FEvalues = join_rows(FEvalues, alphas);
  }
  
  // storage
  List output;
  if (p > 0) {
    output["coefficients"] = coeff;
    output["stderr"] = se;
  }
  output["residuals"] = resid;
  output["FEvalues"] = FEvalues;
  output["alphas"] = computeFixedEffects(data_bak, coeff, FE);
  output["fuck"] = data_bak.col(0) - data_bak.cols(1, p) * coeff;
  if (FEcoefs == 1) {
    //output["mu"] = mu ;
    output["ngroups"] = gtot;
  }
  return(output);
  
}


// [[Rcpp::export()]]
arma::colvec fastplm_predict(double mu,
                             arma::mat FEvalues, // 3 columns
                             arma::mat FE,
                             arma::mat newx,
                             arma::mat beta) {
  // parse data
  int n = newx.n_rows;
  int p = newx.n_cols;
  int m = FE.n_cols; // number of fixed effects
  int g = FEvalues.n_rows;
  
  // number of FE levels in each grouping
  arma::ivec nlvls(m, arma::fill::zeros);
  for (int i=0; i<m; i++) {
    for (int j=0; j<g; j++) {
      if (FEvalues(j,0)==i) {
        nlvls(i) = nlvls(i) + 1;
      }
    }
  }
  
  // starting point of a group of FE coefficients
  arma::ivec gp_start(m, arma::fill::zeros);
  for (int i=1; i<m; i++) {
    gp_start(i) = gp_start(i-1) + nlvls(i-1);
  }
  
  // grand mean
  arma::colvec pred_y(n);
  pred_y.fill(mu);
  
  // covariates
  if (p > 0) {
    pred_y = pred_y + arma::vectorise(newx * beta);
  }
  
  //fixed effects
  for (int i=0; i<n; i++) {
    int check = 0; // check all FE coefficients are found
    for (int j=0; j<m; j++) {
      int cont = 1;
      int lvl = gp_start(j);
      while ((cont == 1) & (lvl<(gp_start(j)+nlvls(j)))) {
        if (FE(i,j) == FEvalues(lvl,1)) {
          pred_y(i) = pred_y(i) + FEvalues(lvl,2);
          cont = 0;
          check++;
        }
        lvl++;
      }
    }
    if (check < m) {
      pred_y(i) = arma::datum::nan;
    }
  }
  return(pred_y);
}
