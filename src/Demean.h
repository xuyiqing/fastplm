#ifndef FASTPLM_DEMEAN_H
#define FASTPLM_DEMEAN_H

// [[Rcpp::plugins(cpp11)]]

#include <algorithm>
#include <vector>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * Mean Indicator
 *
 * A mean indicator is a row vector that captures the operation of computing means.
 * It has three levels, i.e.
 * - an indicator for a specific effect value in some category;
 * - a vector of indicators ("for category"), each for an effect value in some category;
 * - a vector of "for category" ("for datat set"), each for a category in the data set;
 */
struct MeanIndicator {
  typedef std::pair<arma::rowvec, int> ForEffectValue;
  typedef std::vector<ForEffectValue>  ForCategory;
  typedef std::vector<ForCategory>     ForDataSet;
  
  static ForDataSet make(const arma::mat& fixedEffects);
  static void demean(arma::mat& data, const ForCategory& category);
};

enum DemeanTransform {
  kHalperin, kSymmetricHalperin, kCimmino
};

void halperin(arma::mat& data, const MeanIndicator::ForDataSet& categories);
void symmetricHalperin(arma::mat& data, const MeanIndicator::ForDataSet& categories);
void cimmino(arma::mat& data, const MeanIndicator::ForDataSet& categories);

const int kMaxDemeanIteration = 50;
const DemeanTransform kDefaultDemeanTransform = kHalperin;

template<DemeanTransform transform = kDefaultDemeanTransform>
arma::mat demeanMAP(arma::mat data, const arma::mat& fixedEffects) {
  auto meanIndicators = MeanIndicator::make(fixedEffects);
  auto copy = data;
  for (int i = 0; i < kMaxDemeanIteration; i ++) {
    switch (transform) {
    case kHalperin:
      halperin(data, meanIndicators);
      break;
    case kSymmetricHalperin:
      symmetricHalperin(data, meanIndicators);
      break;
    case kCimmino:
      cimmino(data, meanIndicators);
      break;
    default:
      throw std::logic_error("Unknown transform");
    }
    
    auto diff = arma::accu(abs(data - copy));
    if (diff < 1e-5)
      break;
    copy = data;
  }
  return data;
}

#endif
