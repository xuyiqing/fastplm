#include "GPanelMAP.h"

GPanelMAP::GPanelMAP(const arma::cube& X, const arma::mat& Y,
                     const arma::mat& timeInfo, const arma::colvec& indivInfo):
    paramCount(X.n_slices), timeCount(X.n_rows), indivCount(X.n_cols),
    X(X), Y(Y),
    timeInfo(timeInfo), timeParamCount(timeInfo.n_rows),
    indivInfo(indivInfo), indivParamCount(indivInfo.n_rows) {}

void GPanelMAP::buildRegressors() {
    arma::mat augmented;
    arma::rowvec summer;
    
    augmented = timeInfo;
    augmented.insert_cols(0, 1);
    augmented.col(0) = arma::ones(timeCount);
    summer = arma::ones(1 + timeParamCount);
    timeRegressors =
        summer * arma::inv(augmented.t() * augmented) * augmented.t();

    augmented = indivInfo;
    augmented.insert_cols(0, 1);
    augmented.col(0) = arma::ones(indivCount);
    summer = arma::ones(1 + indivParamCount);
    indivRegressors =
        summer * arma::inv(augmented.t() * augmented) * augmented.t();


}

void GPanelMAP::MAP() {
    arma::cube data = arma::join_slices(Y, X);
    arma::cube copy;
    
    do {
        copy = data;
        
        for (int i = 0; i < data.n_slices; i ++) {
            arma::mat panel = data.slice(i);
            
            for (int j = 0; j < indivCount; j ++) {
                arma::colvec col = panel.col(j);
                double delta = arma::as_scalar(timeRegressors * col);
                panel.col(j) -= delta * arma::ones(timeCount);
            }
            
            for (int j = 0; j < timeCount; j ++) {
                arma::colvec col = panel.row(j).t();
                double delta = arma::as_scalar(indivRegressors * col);
                panel.row(j) -= delta * arma::ones(indivCount);
            }
        }
    } while (arma::accu(abs(data - copy)) > 1e-5);
    
    Y_ = data.slice(0);
    X_ = data.slices(1, data.n_slices - 1);
}

arma::colvec GPanelMAP::compute() {
    buildRegressors();
    MAP();
    
    const arma::mat matX(X_.memptr(), timeCount * indivCount, paramCount, false);
    const arma::vec vecY(Y_.memptr(), timeCount * indivCount, 1, false);
    
    return solve(matX, vecY);
}
