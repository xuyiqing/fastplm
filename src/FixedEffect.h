#ifndef FASTPLM_FIXED_EFFECT_H
#define FASTPLM_FIXED_EFFECT_H

#include "Common.h"
#include "Indicator.h"

template <typename Influence>
struct FixedEffect {
    const Indicator& indicator;
    Influence influence;

    typename Influence::SumCollection createSum() const;
    typename Influence::SumItem& get(std::size_t col, typename Influence::SumCollection& sums) const;
    typename Influence::SumItem embed(std::size_t col, double x) const;
    void rescale(typename Influence::SumCollection& sums) const;
    double collapse(std::size_t col, const typename Influence::SumItem& x) const;

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
    typedef arma::vec SumCollection;
};

template<> inline arma::vec FixedEffect<SimpleInfluence>::createSum() const { return arma::zeros(indicator.levelCount); }
template<> inline double& FixedEffect<SimpleInfluence>::get(std::size_t col, arma::vec& sums) const { return sums[col]; }
template<> inline double FixedEffect<SimpleInfluence>::embed(std::size_t col, double x) const { return x; }
template<> inline void FixedEffect<SimpleInfluence>::rescale(arma::vec& sums) const { sums = sums / indicator.levelSizes; }
template<> inline double FixedEffect<SimpleInfluence>::collapse(std::size_t col, const double& x) const { return x; }

typedef FixedEffect<SimpleInfluence> SimpleFixedEffect;

#endif
