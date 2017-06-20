#include <iostream>
#include <fstream>

#include "MatrixReader.h"
#include "GradientDescent.h"
#include "MAP.h"

using namespace std;

int main () {
    auto inputFile = "/Users/selveskii/Documents/Z171.Quantitive Analysis/Projects/fastplm/CXXCode/tests/basic-data.txt";
    ifstream input(inputFile);
    
    auto rawData = readMatrix(input, false);
    arma::mat X = rawData.cols(1, rawData.n_cols - 1);
    arma::mat Y = rawData.col(0);
    
    cout << X << endl;
    cout << Y << endl;
    
    auto rawFixedEffects = readMatrix(input, false);
    std::vector<FixedEffect> fixedEffects;
    for (int i = 0; i < rawFixedEffects.n_cols; i ++)
        fixedEffects.push_back(FixedEffect::fromColumn(rawFixedEffects.col(i)));
    
//    GradientDescent solver1(X, Y, fixedEffects);
    MAP solver2(X, Y, fixedEffects);
//    solver1.compute();
    solver2.compute();
    
    ofstream output("/Users/selveskii/Desktop/result.txt");
//    output << solver1.result.params << endl;
//    output << solver2.result.params << endl;
    cout << solver2.result.params << endl;

//    for (int i = 0; i < solver1.result.effects.size(); i ++)
//        cout << solver1.result.effects[i] << endl;
    
//    for (int i = 0; i < solver2.result.effects.size(); i ++)
//        cout << solver2.result.effects[i] << endl;
    return 0;
}
