#ifndef FASTPLM_GPANEL_MAP_H
#define FASTPLM_GPANEL_MAP_H

#include <vector>

#include "Common.h"

struct GPSolver {
private:
    const arma::uword paramCount;
    const arma::uword timeCount;
    const arma::uword indivCount;
    
    // X[param][time][indiv]
    arma::cube X;
    arma::mat Y;
    
    arma::mat tois;
    arma::mat iots;
    
    const bool isBalanced;
    
    void MAP();
    
    template <bool IsBalanced>
    friend struct BalanceManager;
    
public:
    GPSolver(arma::cube X, arma::mat Y, arma::mat tois, arma::mat iots,
              bool isBalanced = true, bool withFixedEffects = true);
    
    arma::vec compute();
};

#endif
