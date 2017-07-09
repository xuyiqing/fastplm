#include "MatrixReader.h"

template <typename EntryT>
arma::Mat<EntryT> genericReadMatrix(int row, int col, std::istream& input,
                                    bool isColumnMajor) {
    arma::Mat<EntryT> matrix(row, col);
    if (isColumnMajor) {
        for (int i = 0; i < col; i ++)
            for (int j = 0; j < row; j ++)
                input >> matrix(j, i);
    } else {
        for (int i = 0; i < row; i ++)
            for (int j = 0; j < col; j ++)
                input >> matrix(i, j);
    }
    
    return matrix;
}

template <typename EntryT>
arma::Mat<EntryT> genericReadMatrix(std::istream& input, bool isColumnMajor) {
    int row, col;
    input >> row >> col;
    
    return genericReadMatrix<EntryT>(row, col, input, isColumnMajor);
}

arma::mat readMatrix(int row, int col, std::istream& input, bool isColumnMajor) {
    return genericReadMatrix<double>(row, col, input, isColumnMajor);
}

arma::mat readMatrix(std::istream& input, bool isColumnMajor) {
    return genericReadMatrix<double>(input, isColumnMajor);
}
