#ifndef FASTPLM_FIXED_EFFECT_H
#define FASTPLM_FIXED_EFFECT_H

#include "Common.h"
#include "Indicator.h"

template <typename Influence>
struct FixedEffect {
    const Indicator& indicator;
    Influence influence;

    typename Influence::SumCollection createSum() const;
    typename Influence::SumItemRef get(std::size_t col, typename Influence::SumCollection& sums) const;
    typename Influence::SumItem embed(std::size_t col, double x) const;
    void rescale(typename Influence::SumCollection& sums) const;
    double collapse(std::size_t col, typename Influence::SumItemRef x) const;

    typename Influence::SumCollection demean(arma::subview_col<double> data) const {
        typename Influence::SumCollection sums(createSum());

        for (auto i = 0u; i < data.n_rows; i ++)
            get(indicator.indicator[i], sums) += embed(i, data[i]);

        rescale(sums);

        for (auto i = 0u; i < data.n_rows; i ++)
            data[i] -= collapse(i, get(indicator.indicator[i], sums));

        return sums;
    }

    FixedEffect(const Indicator& indicator, Influence&& influence): indicator(indicator), influence(std::move(influence)) {}
};

struct SimpleInfluence {
    typedef double SumItem;
    typedef double& SumItemRef;
    typedef arma::vec SumCollection;
};

typedef FixedEffect<SimpleInfluence> SimpleFixedEffect;

template<> inline arma::vec FixedEffect<SimpleInfluence>::createSum() const {
    return arma::zeros(indicator.levelCount);
}

template<> inline double& FixedEffect<SimpleInfluence>::get(std::size_t col, arma::vec& sums) const {
    return sums[col];
}

template<> inline double FixedEffect<SimpleInfluence>::embed(std::size_t col, double x) const {
    return x;
}

template<> inline void FixedEffect<SimpleInfluence>::rescale(arma::vec& sums) const {
    sums = sums / indicator.levelSizes;
}

template<> inline double FixedEffect<SimpleInfluence>::collapse(std::size_t col, double& x) const {
    return x;
}

struct ComplexInfluence {
    const Indicator& indicator;
    std::size_t dimension;
    arma::mat weight;
    std::vector<arma::mat> normalizers;

    typedef arma::vec SumItem;
    typedef arma::subview_col<double> SumItemRef;
    typedef arma::mat SumCollection;

    ComplexInfluence(const Indicator& indicator, arma::mat&& weight, std::vector<arma::mat>&& normalizers): indicator(indicator), dimension(weight.n_rows), weight(std::move(weight)), normalizers(std::move(normalizers)) {}
};

typedef FixedEffect<ComplexInfluence> ComplexFixedEffect;

template<> inline arma::mat FixedEffect<ComplexInfluence>::createSum() const {
    return arma::zeros(influence.dimension, indicator.levelCount);
}

template<> inline arma::subview_col<double> FixedEffect<ComplexInfluence>::get(std::size_t i, arma::mat& sums) const {
    return sums.col(i);
}

template<> inline arma::vec FixedEffect<ComplexInfluence>::embed(std::size_t col, double x) const {
    const auto level = influence.indicator.indicator[col];
    return influence.weight.col(level) * x;
}

template<> inline void FixedEffect<ComplexInfluence>::rescale(arma::mat& sums) const {
    for (auto i = 0; i < sums.n_cols; i ++)
        sums.col(i) = influence.normalizers[i] * sums.col(i);
}

template<> inline double FixedEffect<ComplexInfluence>::collapse(std::size_t col, arma::subview_col<double> x) const {
    const auto level = influence.indicator.indicator[col];
    return arma::dot(x, influence.weight.col(level));
}

#endif
