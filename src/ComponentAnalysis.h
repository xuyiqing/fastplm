#ifndef FASTPLM_COMPONENT_ANALYSIS_H
#define FASTPLM_COMPONENT_ANALYSIS_H

#include <unordered_set>

#include "Common.h"
#include "Indicator.h"

typedef std::vector<arma::uvec> ComponentTables;

optional<ComponentTables> computeComponents(const arma::uvec& groupSizes, const std::vector<Indicator>& indicators);

struct CrossComponentError {
    std::size_t groupX, groupY;
    std::size_t valueX, valueY;
    std::size_t row;

#ifndef BUILD_WITHOUT_R
    operator Rcpp::List() const {
        Rcpp::List _;
        _["group.x"] = groupX + 1;
        _["group.y"] = groupY + 1;

        _["value.x"] = valueX + 1;
        _["value.y"] = valueY + 1;

        _["row"] = row + 1;

        return _;
    }
#endif

};

std::vector<CrossComponentError> checkComponents(const ComponentTables& tables, const arma::umat& indicators);

#endif
