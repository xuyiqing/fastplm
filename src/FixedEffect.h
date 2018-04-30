#ifndef FASTPLM_FIXED_EFFECT_H
#define FASTPLM_FIXED_EFFECT_H

#include "Common.h"
#include "Indicator.h"

template <typename InfluenceType, typename SumItem, typename SumCollection>
struct FixedEffect {
    const Indicator& indicator;
    InfluenceType influence;

    SumCollection createSum() const;
    SumItem& get(std::size_t col, SumCollection& sums) const;
    SumItem embed(std::size_t col, double x) const;
    void rescale(SumCollection& sums) const;
    double collapse(std::size_t col, const SumItem& x) const;

    SumCollection demean(arma::subview_col<double> data) const {
        SumCollection sums(createSum());

        for (auto i = 0u; i < data.n_rows; i ++)
            get(indicator.indicator[i], sums) += embed(i, data[i]);

        rescale(sums);

        for (auto i = 0u; i < data.n_rows; i ++)
            data[i] -= collapse(i, get(indicator.indicator[i], sums));

        return sums;
    }

    FixedEffect(const Indicator& indicator, InfluenceType&& influence): indicator(indicator), influence(std::move(influence)) {}
};

enum Singleton { singleton };

template<> inline arma::vec FixedEffect<Singleton, double, arma::vec>::createSum() const { return arma::zeros(indicator.levelCount); }
template<> inline double& FixedEffect<Singleton, double, arma::vec>::get(std::size_t col, arma::vec& sums) const { return sums[col]; }
template<> inline double FixedEffect<Singleton, double, arma::vec>::embed(std::size_t col, double x) const { return x; }
template<> inline void FixedEffect<Singleton, double, arma::vec>::rescale(arma::vec& sums) const { sums = sums / indicator.levelSizes; }
template<> inline double FixedEffect<Singleton, double, arma::vec>::collapse(std::size_t col, const double& x) const { return x; }

typedef FixedEffect<Singleton, double, arma::vec> SimpleFixedEffect;

#endif
