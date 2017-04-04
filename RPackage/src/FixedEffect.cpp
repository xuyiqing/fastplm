#include "FixedEffect.h"

FixedEffect::FixedEffect(const arma::mat& indicator_, FixedEffect::Indices indices)
: indicator(indicator_), groupCount(indicator_.n_cols), observationCount(indicator_.n_rows), indices(indices)
{
    valuesOccurences = std::vector<int>(groupCount);
    for (arma::uword i = 0; i < groupCount; i ++) {
        int count = 0;
        auto col = indicator.col(i);
        for (arma::uword j = 0; j < observationCount; j ++)
            if (col(j) != 0)
                count ++;
        valuesOccurences[i] = count;
    }
}

FixedEffect FixedEffect::fromColumn(const arma::colvec& column, FixedEffect::Indices indices) {
    size_t groupCount = indices.size();
    arma::mat indicator = arma::zeros(column.n_rows, groupCount);
    
    for (int i = 0; i < column.n_rows; i ++) {
        double effectValue = column(i);
        indicator(i, indices[effectValue]) = 1;
    }
    
    return FixedEffect(indicator, indices);
}

FixedEffect FixedEffect::fromColumn(const arma::colvec& column) {
    arma::colvec effectValues = arma::unique(column);
    std::sort(effectValues.begin(), effectValues.end());
    FixedEffect::Indices indices;
    for (int i = 0; i < effectValues.size(); i ++)
        indices.insert({effectValues(i), i});
    return FixedEffect::fromColumn(column, indices);
}
