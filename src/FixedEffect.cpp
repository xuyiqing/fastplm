#include <unordered_set>
#include <queue>

#include "CrashQueue.h"
#include "FixedEffect.h"

FixedEffect::FixedEffect() {}

FixedEffect FixedEffect::fromColumn(const arma::colvec& column) {
    // We can assume that the column consists of consecutive integers from 1 to N.
    int groupCount = 0;
    
    for (double val : column)
        groupCount = std::max(groupCount, static_cast<int>(val));
    
    FixedEffect effect;
    effect.groupCount = groupCount;
    
    effect.column = std::vector<int>();
    effect.column.reserve(column.n_rows);
    for (int i = 0; i < column.n_rows; i ++)
        effect.column.push_back(static_cast<int>(column[i]) - 1);
    
    effect.groupSizes = std::vector<int>(groupCount, 0);
    for (int val : effect.column)
        effect.groupSizes[val] ++;

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
