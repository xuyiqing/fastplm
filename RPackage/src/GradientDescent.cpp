#include "GradientDescent.h"

GradientDescent::GradientDescent(const arma::mat& X, const arma::colvec& Y,
                                 const std::vector<FixedEffect>& fixedEffects,
                                 bool computeXtXInverse):
    PlainModel(X, Y, fixedEffects, computeXtXInverse) {}

arma::colvec GradientDescent::nextParams(const arma::colvec& partialResidual) {
    return XtXInverse * (X.t() * partialResidual);
}

arma::colvec GradientDescent::nextFixedEffect(const FixedEffect& effect, const arma::colvec& partialResidual) {
    auto groupCount = effect.groupCount;

    arma::colvec next(groupCount);
    for (int i = 0; i < groupCount; i ++) {
        auto sum = arma::dot(effect.indicator.col(i), partialResidual);
        next(i) = sum / effect.valuesOccurences[i];
    }
    next -= (arma::sum(next) / next.n_rows) * arma::ones(next.n_rows);
    return next;
}

double GradientDescent::nextIntercept(const arma::colvec& partialResidual) {
    return arma::sum(partialResidual) / observationCount;
}

arma::colvec GradientDescent::computeResidual(const GradientDescent::Result& result) {
    arma::colvec residual = Y;
    residual -= X * result.params;
    for (int i = 0; i < fixedEffects.size(); i ++) {
        residual -= fixedEffects[i].indicator * result.effects[i];
    }
    residual -= result.intercept * arma::ones(observationCount);
    return residual;
}

void GradientDescent::compute() {
    auto& current = result;
    current.params = arma::zeros(paramCount);
    current.effects = std::vector<arma::colvec>(fixedEffects.size());
    for (int i = 0; i < current.effects.size(); i ++) {
        current.effects[i] = arma::zeros(fixedEffects[i].groupCount);
    }
    current.intercept = 0.0;
    
    Result last;
    do {
        last = current;
        arma::colvec residual = computeResidual(last);
        current.params = nextParams(residual + X * last.params);
        for (int i = 0; i < fixedEffects.size(); i ++) {
            auto partialResidual = residual + fixedEffects[i].indicator * last.effects[i];
            current.effects[i] = nextFixedEffect(fixedEffects[i], partialResidual);
        }
        current.intercept = nextIntercept(residual + last.intercept * arma::ones(observationCount));
    } while (current.computeDiff(last) > 1e-5);
}

double GradientDescent::Result::computeDiff(const GradientDescent::Result& last) {
    double diff = arma::accu(abs(this->params - last.params));
    
    for (int i = 0; i < this->effects.size(); i ++) {
        const double effectDiff = arma::accu(arma::abs(this->effects[i] - last.effects[i]));
        if (effectDiff > diff)
            diff = effectDiff;
    }
    const double interceptDiff = std::abs(this->intercept - last.intercept);
    if (interceptDiff > diff)
        diff = interceptDiff;
    return diff;
}
