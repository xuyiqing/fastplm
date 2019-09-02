#include "ComponentAnalysis.h"

// Below is an implementation of the union-find algorithm with path compression and union by rank.
typedef arma::uword ID;

struct Component {
    ID parent, rank;
};

ID getComponent(std::vector<Component>& components, ID id) {
    if (components[id].parent != id)
        components[id].parent = getComponent(components, components[id].parent);
    return components[id].parent;
}

void unionByRank(std::vector<Component>& components, ID x, ID y) {
    auto cx = getComponent(components, x);
    auto cy = getComponent(components, y);

    if (components[cx].rank < components[cy].rank)
        components[cx].parent = cy;
    else if (components[cx].rank > components[cy].rank)
        components[cy].parent = cx;
    else {
        components[cy].parent = cx;
        components[cx].rank ++;
    }
}
// Above is an implementation of the union-find algorithm with path compression and union by rank.

optional<ComponentTables> computeComponents(const arma::uvec& groupSizes, const std::vector<Indicator>& indicators) {
    auto totalNodes = arma::sum(groupSizes);
    std::vector<Component> components(totalNodes);
    for (auto i = 0u; i < totalNodes; i ++) {
        components[i].parent = i;
        components[i].rank = 0;
    }

    std::size_t effectCount = groupSizes.n_elem;
    std::size_t totalRows = indicators[0].indicator.n_rows;
    for (auto row = 0u; row < totalRows; row ++) {
        ID offset = 0;
        for (auto i = 0u; i < effectCount - 1; i ++) {
            ID x = offset + indicators[i].indicator[row];
            offset += groupSizes[i];
            ID y = offset + indicators[i + 1].indicator[row];
            unionByRank(components, x, y);
        }
    }

    std::unordered_set<ID> distinctGroups;
    ComponentTables tables;
    ID offset = 0;

    std::transform(groupSizes.begin(), groupSizes.end(), std::back_inserter(tables), [&](auto groupSize) {
        arma::uvec table(groupSize);

        for (auto i = 0u; i < groupSize; i ++) {
            auto group = getComponent(components, i + offset);
            distinctGroups.insert(group);
            table[i] = group;
        }

        offset += groupSize;
        return table;
    });

    if (distinctGroups.size() <= 1)
        return std::experimental::nullopt;

    return std::experimental::make_optional(tables);
}

std::vector<CrossComponentError> checkComponents(const ComponentTables& tables, const arma::umat& indicators) {
    std::vector<CrossComponentError> errors;

    for (auto row = 0u; row < indicators.n_rows; row ++) {
        for (auto col = 0u; col < indicators.n_cols - 1; col ++) {
            if (std::isnan(indicators(row, col)))
                break;

            ID x = indicators(row, col);
            ID y = indicators(row, col + 1);

            if (tables[col][x] == tables[col + 1][y])
                continue;

            CrossComponentError error;
            error.groupX = col;
            error.groupY = col + 1;
            error.valueX = x;
            error.valueY = y;
            error.row = row;
            errors.push_back(error);
        }
    }

    return errors;
}
