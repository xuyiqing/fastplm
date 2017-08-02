#include <unordered_set>
#include <queue>

#include "CrashQueue.h"
#include "FixedEffect.h"

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

struct DemeanArgs {
    double* vector;
    const std::vector<int>& ids;
    const std::vector<int>& groupSizes;

    DemeanArgs(double* vector, const std::vector<int>& ids, const std::vector<int>& groupSizes)
    :vector(vector), ids(ids), groupSizes(groupSizes) {}
};

void demeanColumn(void* ptr) {
    DemeanArgs* casted = static_cast<DemeanArgs*>(ptr);
    auto vector = casted->vector;
    auto ids = casted->ids;
    const auto& groupSizes = casted->groupSizes;
    std::vector<double> count(groupSizes.size(), 0.0);

    for (int i = 0; i < ids.size(); i ++)
        count[ids[i]] += vector[i];
    for (int i = 0; i < groupSizes.size(); i ++)
        count[i] /= static_cast<double>(groupSizes[i]);
    for (int i = 0; i < ids.size(); i ++)
        vector[i] -= count[ids[i]];
}

void FixedEffect::demean(arma::mat& data) const {
    std::vector<DemeanArgs> payloads;
    payloads.reserve(data.n_cols);
    for (std::size_t i = 0; i < data.n_cols; i ++)
        payloads.emplace_back(data.colptr(i), column, groupSizes);
    std::queue<void*> casted;
    for (std::size_t i = 0; i < data.n_cols; i ++)
        casted.push(&payloads[i]);

    if (mainQueue == nullptr) {
        for (std::size_t i = 0; i < data.n_cols; i ++)
            demeanColumn(&payloads[i]);
    } else {
        mainQueue->commit(demeanColumn, std::move(casted));
        mainQueue->crash();
    }
}
