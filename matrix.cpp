#include <vector>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include "matrix.h"

namespace matrix {

namespace {

const float FLOAT_GAP = 1e-6f;
const double DOUBLE_GAP = 1e-12;

bool is_match(float a, float b) { //instead of ==
    return (fabsf(a - b) < FLOAT_GAP);
}

bool is_match(double a, double b) { //instead of ==
    return (fabs(a - b) < DOUBLE_GAP);
}

bool is_match(int a, int b) { //instead of ==
    return (a == b);
}

bool is_match(long long int a, int b) { //instead of ==
    return (a == b);
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

//it is support function only for determinant()
template <typename T, typename T_i> void L_step(const Square_Matrix<T>& matrix,
                                                std::vector<std::vector<T_i>>& L,
                                                std::vector<std::vector<T_i>>& U,
                                                std::vector<T_i>& main_diag,
                                                size_t i) {
    L[i - 1][0] = static_cast<T_i>(matrix(i, 0)) / main_diag[0];

    for(size_t j = 1; j < i; ++j) {
        L[i - 1][j] = (static_cast<T_i>(matrix(i, j)) - sum_for_determinant<T_i>(L[i - 1], U[j - 1], j)) / main_diag[j];
    }
}

//it is support function only for determinant()
template <typename T, typename T_i> void U_step(const Square_Matrix<T>& matrix,
                                                std::vector<std::vector<T_i>>& L,
                                                std::vector<std::vector<T_i>>& U,
                                                size_t i) {

    for(size_t j = i; j < (matrix.ret_size() - 1); ++j) {
        U[j][i] = static_cast<T_i>(matrix(i, j + 1)) - sum_for_determinant<T_i>(L[i - 1], U[j], i);
    }
}

template <typename T, typename T_i> T_i universal_determinant(const Square_Matrix<T>& input_matrix) {
    if(input_matrix.ret_column_size() == 0) {
        return static_cast<T_i>(0);
    }

    Square_Matrix<T> matrix = input_matrix;
    const size_t matr_size = matrix.ret_size();
    std::vector<std::vector<T_i>> L;
    std::vector<std::vector<T_i>> U; //for LU decomposition it is transposed U (for speed)
    std::vector<T_i> main_diag; //of U

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


    //It isn't standart LU decompose method, because it should work correctly with some main minor = 0
    //
    //L and U are counted line by line and synchronously by steps
    //At each step we count (in this order):
    //1) 1 ROW of L
    //2) 1 diagonal elem
    //3) CHECK
    //4) 1 COLUMN of TRANSPOSED U (1 ROW of U)
    //
    //CHECK: if we have 0 in U main diag, there are two cases:
    //
    //1. we can fix it by adding other matrix rows under current row
    //   trying to make main minor nondegen
    //
    //   In this case we sholud repeat points 1), 2) and 3) after each attempt
    //
    //2. if (1) way can't fix main minor so the matrix is degen
    //   determinant = 0
    //
    //First step should be different with others (for convenience)
    //point 1) isn't exists there (but it isn't all difference)
    size_t str_n = 1;
    while(is_match(matrix(0, 0), static_cast<T>(0))) {
        if(str_n == matr_size) {
            return static_cast<T_i>(0);
        }
        matrix.add_row_to_row(str_n, 0);
        ++str_n;
    }

    main_diag[0] = static_cast<T_i>(matrix(0, 0));
    for(size_t j = 0; j < matr_size - 1; ++j) {
        U[j][0] = static_cast<T_i>(matrix(0, j + 1));
    }


    //other steps after the first
    for(size_t i = 1; i < matr_size; ++i) {
        size_t str_n;

        //building L
        L_step<T, T_i>(matrix, L, U, main_diag, i);

        //building main_diag, checking two special cases
        main_diag[i] = static_cast<T_i>(matrix(i, i)) - sum_for_determinant<T_i>(L[i - 1], U[i - 1], i);

        //check
        str_n = i + 1;
        while(is_match(main_diag[i], static_cast<T_i>(0))) {
            if(str_n == matr_size) { //unfixing case
                return static_cast<T_i>(0);
            }

            //fix matrix
            matrix.add_row_to_row(str_n, i);

            //rebuild L row and main_diag elem
            L_step<T, T_i>(matrix, L, U, main_diag, i);
            main_diag[i] = static_cast<T_i>(matrix(i, i)) - sum_for_determinant<T_i>(L[i - 1], U[i - 1], i);

            ++str_n;
        }


        //building U
        U_step<T, T_i>(matrix, L, U, i);

    }

    //calculating determinant
    T_i det = main_diag[0];
    for(size_t i = 1; i < matr_size; ++i) {
        det *= main_diag[i];
    }
    return det;
}

} //ending namespace {}

double determinant(const Square_Matrix<double>& input_matrix) {
    return universal_determinant<double, double>(input_matrix);
}

float determinant(const Square_Matrix<float>& input_matrix) {
    return static_cast<float>(universal_determinant<float, double>(input_matrix));
}

int determinant(const Square_Matrix<int>& input_matrix) {
    return lround(universal_determinant<int, double>(input_matrix));
}

long long int determinant(const Square_Matrix<long long int>& input_matrix) {
    return llround(universal_determinant<long long int, double>(input_matrix));
}

//... other types if will be need

} //ending namespace matrix {}
