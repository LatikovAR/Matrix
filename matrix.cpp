#include <vector>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include "matrix.h"

namespace matrix {

namespace {

const double DOUBLE_GAP = 0.000001;

bool is_match(double a, double b) { //== for double
    return (fabs(a - b) < DOUBLE_GAP);
}

//it is support function only for determinant()
template <typename T> T sum_for_determinant(std::vector<T>& L_str,
                                            std::vector<T>& U_str,
                                            size_t n) {
    assert(n <= L_str.size());
    assert(n <= U_str.size());
    assert(n >= 0);

    T sum = static_cast<T>(0);
    for(size_t i = 0; i < n; ++i) {
        sum += L_str[i] * U_str[i];
    }
    return sum;
}

//this function can change matrix
//external functions should give a copy in it
template <typename T> T universal_determinant(Square_Matrix<T>& matrix) {
    if(matrix.ret_column_size() == 0) {
        return static_cast<T>(0);
    }

    const size_t matr_size = matrix.ret_column_size();
    std::vector<std::vector<T>> L;
    std::vector<std::vector<T>> U; //for LU decomposition it is transposed U (for speed)
    std::vector<T> main_diag; //of U

    //L and transposed U store only elems under their main diagonals
    //main_diag - it's U main diagonal
    //L main diagonal should have only "1" elements and not stores
    L.resize(matr_size - 1);
    for(size_t i = 0; i < (matr_size - 1); ++i) {
        L[i].resize(i + 1);
    }
    U.resize(matr_size - 1);
    for(size_t i = 0; i < (matr_size - 1); ++i) {
        U[i].resize(i + 1);
    }
    main_diag.resize(matr_size);


    //first step should be different with others
    //
    //if we have 0 in U main diag, there are two cases:
    //
    //1. we can fix it by adding other matrix rows under current row
    //   trying to make main minor nondegen
    //
    //2. if (1) way can't fix main minor so the matrix is degen
    //   determinant = 0
    //
    //these two cases also checked not only at the first step
    size_t str_n = 1;
    while(is_match(matrix(0, 0), 0)) { //TO DO: operator==
        if(str_n == matr_size) {
            return static_cast<T>(0);
        }
        matrix.add_row_to_row(str_n, 0);
        ++str_n;
    }

    main_diag[0] = matrix(0, 0);
    for(size_t j = 0; j < matr_size - 1; ++j) {
        U[j][0] = matrix(0, j + 1);
    }
    for(size_t j = 0; j < matr_size - 1; ++j) {
        L[j][0] = matrix(j + 1, 0) / main_diag[0];
    }


    //other steps after the first
    for(size_t i = 1; i < matr_size; ++i) {
        size_t str_n;
        T buf;

        //building main_diag, checking two special cases
        buf = sum_for_determinant<T>(L[i - 1], U[i - 1], i);
        str_n = i + 1;
        main_diag[i] = matrix(i, i) - buf;
        while(is_match(main_diag[i], 0)) {
            if(str_n == matr_size) {
                return static_cast<T>(0);
            }
            matrix.add_row_to_row(str_n, i);
            main_diag[i] = matrix(i, i) - buf;
            ++str_n;
        }


        //building U
        for(size_t j = i; j < (matr_size - 1); ++j) {
            U[j][i] = matrix(i, j + 1) - sum_for_determinant<T>(L[i - 1], U[j], i);
        }

        //building L
        for(size_t j = i; j < (matr_size - 1); ++j) {
            L[j][i] = (matrix(j + 1, i) - sum_for_determinant<T>(L[j], U[i - 1], i)) / main_diag[i];
        }
    }

    //calculating determinant
    T det = main_diag[0];
    for(size_t i = 1; i < matr_size; ++i) {
        det *= main_diag[i];
    }
    return det;
}

} //ending namespace {}

double determinant(const Square_Matrix<double>& input_matrix) {
    Square_Matrix<double> matrix = input_matrix;
    double det = universal_determinant(matrix);
    return det;
}

int determinant(const Square_Matrix<int>& input_matrix) {
    //Square_Matrix<int> -> Square_Matrix<double>
    int det = static_cast<int>(universal_determinant(matrix));
    return det;
}

} //ending namespace matrix {}