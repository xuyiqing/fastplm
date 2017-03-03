#include "Demean.h"

// MeanIndicator implementations.
MeanIndicator::ForDataSet MeanIndicator::make(const arma::mat& fixedEffects) {
  MeanIndicator::ForDataSet forDataSet;
  fixedEffects.each_col([&](const arma::colvec fixedEffect) {
    arma::colvec effectValues = arma::unique(fixedEffect);
    MeanIndicator::ForCategory forCategory;
    effectValues.for_each([&](double effectValue) {
      arma::rowvec indicator = arma::zeros(1, fixedEffect.n_rows);
      int count = 0;
      for (int i = 0; i < fixedEffect.n_rows; i ++) {
        if (fixedEffect(i) == effectValue) {
          indicator(i) = 1;
          count ++;
        }
      }
      forCategory.push_back(std::make_pair(indicator, count));
    });
    forDataSet.push_back(forCategory);
  });
  return forDataSet;
}

void MeanIndicator::demean(arma::mat& data, const MeanIndicator::ForCategory& category) {
  for (const auto& forEffectValue : category) {
    const auto& indicator = forEffectValue.first;
    int count = forEffectValue.second;
    if (count > 0) {
      arma::rowvec means = (indicator * data) / static_cast<double>(count);
      data -= (indicator.t() * means);
    }
  }
}

// DemeanTransformation implementations.
void halperin(arma::mat& data, const MeanIndicator::ForDataSet& categories) {
  for (const auto& category : categories)
    MeanIndicator::demean(data, category);
}

void symmetricHalperin(arma::mat& data, const MeanIndicator::ForDataSet& categories) {
  int i = 0;
  for (; i < categories.size() - 1; i ++)
    MeanIndicator::demean(data, categories[i]);
  for (; i >= 0; i --)
    MeanIndicator::demean(data, categories[i]);
}

void cimmino(arma::mat& data, const MeanIndicator::ForDataSet& categories) {
  arma::mat overallMeans = arma::zeros(data.n_rows, data.n_cols);
  
  for (const auto& category : categories)
    for (const auto& forEffectValue : category) {
      const auto& indicator = forEffectValue.first;
      int count = forEffectValue.second;
      if (count > 0) {
        arma::rowvec means = (indicator * data) / static_cast<double>(count);
        overallMeans += (indicator.t() * means);
      }
    }
    data -= overallMeans / categories.size();
}
