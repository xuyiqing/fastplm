#include "MatrixReader.h"

template <typename EntryT>
arma::Mat<EntryT> genericReadMatrix(std::istream& input) {
    int row, column;
    input >> row >> column;
    
    arma::Mat<EntryT> matrix(row, column);
    for (int i = 0; i < column; i ++)
        for (int j = 0; j < row; j ++)
            input >> matrix(j, i);
    return matrix;
}

arma::mat readMatrix(std::istream& input) {
    return genericReadMatrix<double>(input);
}
