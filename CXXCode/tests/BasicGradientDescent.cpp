#include <iostream>
#include <fstream>

#include "MatrixReader.h"
#include "GradientDescent.h"

using namespace std;

int main () {
    auto inputFile = "/Users/selveskii/Documents/Z171.Quantitive Analysis/Projects/better-fastplm/CXXCode/tests/basic-data.txt";
    ifstream input(inputFile);
    
    auto rawData = readMatrix(input);
    arma::mat X = rawData.cols(1, rawData.n_cols - 1);
    arma::mat Y = rawData.col(0);
    
    auto rawFixedEffects = readMatrix(input);
    std::vector<FixedEffect> fixedEffects;
    for (int i = 0; i < rawFixedEffects.n_cols; i ++)
        fixedEffects.push_back(FixedEffect::fromColumn(rawFixedEffects.col(i)));
    
    GradientDescent solver(X, Y, fixedEffects);
    solver.compute();
    auto result = solver.result;
    
    ofstream output("/Users/selveskii/Desktop/result.txt");
    output << result.params << endl;
    output << result.intercept << endl;
    return 0;
}
