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
    effect.observationCount = column.n_rows;
    
    effect.valuesOccurences = std::vector<int>(effect.groupCount);
    std::transform(listOfValueCounters.cbegin(), listOfValueCounters.cend(),
                   effect.valuesOccurences.begin(),
                   [](std::pair<double, int> x) { return x.second; } );
    
    effect.indices = FixedEffect::Indices();
    for (int i = 0; i < listOfValueCounters.size(); i ++)
        effect.indices.insert({listOfValueCounters[i].first, i});
    
    effect.indicator = arma::zeros(column.n_rows, effect.groupCount);
    for (int i = 0; i < column.n_rows; i ++) {
        double effectValue = column(i);
        effect.indicator(i, effect.indices[effectValue]) = 1;
    }
    
    return effect;
}
