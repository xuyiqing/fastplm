#ifndef FASTPLM_COMMON_H
#define FASTPLM_COMMON_H

#include <array>
#include <vector>

#ifdef BUILD_WITHOUT_R

#include <armadillo>

#elif

// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

#endif
#endif
