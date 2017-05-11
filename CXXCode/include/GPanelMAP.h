#ifndef FASTPLM_GPANEL_MAP_H
#define FASTPLM_GPANEL_MAP_H

#include "Common.h"

struct GPanelMAP {
private:
    const arma::uword paramCount;
    const arma::uword timeCount;
    const arma::uword indivCount;
    
    // X[param][time][indiv]
    const arma::cube X;
    const arma::mat Y;
    
    // timeInfo[time][param]
    const arma::mat timeInfo;
    const arma::uword timeParamCount;
    
    const arma::mat indivInfo;
    const arma::uword indivParamCount;
    
public:
    GPanelMAP(const arma::cube& X, const arma::mat& Y,
              const arma::mat& timeInfo, const arma::colvec& indivInfo);

    arma::colvec compute();
    
private:
    arma::mat timeRegressors;
    arma::mat indivRegressors;
    
    arma::cube X_;
    arma::mat Y_;
    
    void buildRegressors();
    void MAP();
};

#endif
