#include "CrushQueue.h"
#include "FixedEffects.h"
#include "FixedEffectModel.h"

int main () {
    arma::mat data;
    data.load("data.csv", arma::file_type::csv_ascii);
    arma::mat inds;
    inds.load("inds.csv", arma::file_type::csv_ascii);

    arma::uvec levelCounts = { 50, 50, 50 };
    arma::uvec simpleEffects = { 0, 1, 2 };
    arma::uvec zeroVec;

    auto fes = FixedEffects::create(levelCounts, inds, simpleEffects, zeroVec, zeroVec, std::vector<arma::mat>());

    mainQueue = new CrushQueue(1);
    ScopeGuard _([]{ delete mainQueue; mainQueue = nullptr; });

    auto model = FixedEffectModel::solve(data, *fes);

    std::cout << model.demeaned.beta << std::endl;

    return 0;
}
