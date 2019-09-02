#ifndef FASTPLM_COMMON_H
#define FASTPLM_COMMON_H

#include <array>
#include <vector>
#include <functional>
#include <experimental/optional>

using std::experimental::optional;



class ScopeGuard {
private:
    std::function<void()> toCall;
public:
    ScopeGuard(std::function<void()> toCall): toCall(toCall) {}
    ~ScopeGuard() {
        toCall();
    }
};

const double DOUBLE_TOLERANCE = 1e-5;

#ifdef BUILD_WITHOUT_R

#include <armadillo>

#else

// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#endif

inline bool isZero(const arma::colvec& vec, const double eps = DOUBLE_TOLERANCE) {
    return arma::all(arma::abs(vec) < eps);
}

inline const arma::colvec getY(const arma::mat& data) {
    return data.col(0);
}

inline const arma::mat getX(const arma::mat& data) {
    if (data.n_cols > 1)
        return data.cols(1, data.n_cols - 1);
    else
        return arma::mat( data.n_rows, 0 );
}

#endif
