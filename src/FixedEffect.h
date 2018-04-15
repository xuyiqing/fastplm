#ifndef FASTPLM_FIXED_EFFECT_H
#define FASTPLM_FIXED_EFFECT_H

#include "Common.h"

template <typename InfluenceType, typename SumItem, typename SumCollection>
struct FixedEffect {
    std::size_t levelCount;
    const arma::subview_col<arma::uword> indicators;
    InfluenceType influence;

    SumCollection createSum() const;
    SumItem& get(std::size_t col, SumCollection& sums) const;
    SumItem embed(std::size_t col, double x) const;
    void rescale(SumCollection& sums) const;
    double collapse(std::size_t col, const SumItem& x) const;

    SumCollection demean(arma::subview_col<double> data) const {
        SumCollection sums(createSum());

        for (auto i = 0u; i < data.n_rows; i ++)
            get(indicators[i], sums) += embed(i, data[i]);

        rescale(sums);

        for (auto i = 0u; i < data.n_rows; i ++)
            data[i] -= collapse(i, get(indicators[i], sums));

        return sums;
    }

    FixedEffect(std::size_t levelCount, const arma::subview_col<arma::uword> indicators, InfluenceType&& influence): levelCount(levelCount), indicators(indicators), influence(std::move(influence)) {}
};

template<> inline arma::vec FixedEffect<arma::vec, double, arma::vec>::createSum() const { return arma::zeros(levelCount); }
template<> inline double& FixedEffect<arma::vec, double, arma::vec>::get(std::size_t col, arma::vec& sums) const { return sums[col]; }
template<> inline double FixedEffect<arma::vec, double, arma::vec>::embed(std::size_t col, double x) const { return x; }
template<> inline void FixedEffect<arma::vec, double, arma::vec>::rescale(arma::vec& sums) const { sums = sums / influence; }
template<> inline double FixedEffect<arma::vec, double, arma::vec>::collapse(std::size_t col, const double& x) const { return x; }

typedef FixedEffect<arma::vec, double, arma::vec> SimpleFixedEffect;

#endif
