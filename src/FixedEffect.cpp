#include "FixedEffect.h"

#include <unordered_set>

FixedEffect::FixedEffect() {}

FixedEffect FixedEffect::fromColumn(const arma::colvec& column) {
    std::unordered_map<double, int> valueCounters;
    
    for (double val : column) {
        auto iter = valueCounters.find(val);
        if (iter == valueCounters.end())
            valueCounters.insert({ val, 1 });
        else
            iter->second ++;
    }
    
    std::vector<std::pair<double, int>>
    listOfValueCounters(valueCounters.begin(), valueCounters.end());
    std::sort(listOfValueCounters.begin(), listOfValueCounters.end(),
              [](std::pair<double, int> lhs, std::pair<double, int> rhs){
                  return lhs.first < rhs.first;
              });
    
    FixedEffect effect;
    effect.groupCount = valueCounters.size();
    
    effect.groupSizes = std::vector<int>(effect.groupCount);
    std::transform(listOfValueCounters.cbegin(), listOfValueCounters.cend(),
                   effect.groupSizes.begin(),
                   [](std::pair<double, int> x) { return x.second; } );
    
    effect.indices = std::unordered_map<double, int>();
    for (int i = 0; i < listOfValueCounters.size(); i ++)
        effect.indices.insert({listOfValueCounters[i].first, i});
    
    effect.column = std::vector<int>();
    effect.column.reserve(column.n_rows);
    for (int i = 0; i < column.n_rows; i ++) {
        double effectValue = column(i);
        effect.column.push_back(effect.indices[effectValue]);
    }
    
    return effect;
}

std::vector<double> FixedEffect::computeMean(const double* ptr) const {
    std::vector<double> count(groupCount, 0.0);
    
    for (int i = 0; i < column.size(); i ++)
        count[column[i]] += ptr[i];
    for (int i = 0; i < groupCount; i ++)
        count[i] /= static_cast<double>(groupSizes[i]);
    
    return count;
}

void FixedEffect::demean(arma::mat& data) const {
    for (int i = 0; i < data.n_cols; i ++) {
        std::vector<double> count(groupCount, 0.0);
        
        double *ptr = data.colptr(i);
        for (int i = 0; i < column.size(); i ++)
            count[column[i]] += ptr[i];
        for (int i = 0; i < groupCount; i ++)
            count[i] /= static_cast<double>(groupSizes[i]);
        
        for (int i = 0; i < column.size(); i ++)
            ptr[i] -= count[column[i]];
    }
}
