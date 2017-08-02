#ifndef FASTPLM_COMMON_H
#define FASTPLM_COMMON_H

#include <array>
#include <vector>
#include <functional>

class ScopeGuard {
private:
    std::function<void()> toCall;
public:
    ScopeGuard(std::function<void()> toCall): toCall(toCall) {}
    ~ScopeGuard() {
        toCall();
    }
};

#ifdef BUILD_WITHOUT_R

#include <armadillo>

#else

// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

#endif
#endif
