#ifndef FASTPLM_MATRIX_READER_H
#define FASTPLM_MATRIX_READER_H

#include <istream>

#include "Common.h"

arma::mat readMatrix(int row, int col, std::istream& input,
                     bool isColumnMajor = true);
arma::mat readMatrix(std::istream& input,
                     bool isColumnMajor = true);

#endif
