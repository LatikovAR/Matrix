#include <iostream>
#include <cstdlib>
#include <cassert>
#include <ctime>

#include "matrix.h"

using namespace matrix;

void test0() {
    Square_Matrix<int> m({{1, 1, 1},
                          {1, 1, 1},
                          {1, 1, 1}});

    int det = determinant(m);

    assert(det == 0);
    std::cout << "Test0 complete.\n";
}

void test1() {
    Square_Matrix<int> m({{2, 0, 1},
                          {1, 0, 1},
                          {1, 1, 1}});

    int det = determinant(m);

    assert(det == -1);
    std::cout << "Test1 complete.\n";
}

void test2() {
    Square_Matrix<int> m({{0, 0, 1, 1},
                          {0, 0, 1, 2},
                          {1, 1, 1, 1},
                          {1, 2, 1, 1}});

    int det = determinant(m);

    assert(det == 1);
    std::cout << "Test2 complete.\n";
}

namespace  {

Square_Matrix<long long int> make_add_matrix(const Square_Matrix<long long int>& matrix,
                                             size_t e_row,
                                             size_t e_column) {
    assert(e_row < matrix.size());
    assert(e_column < matrix.size());
    std::vector<std::vector<long long int>> new_rows;
    new_rows.resize(matrix.size() - 1);
    for(size_t i = 0; i < (matrix.size() - 1); ++i) {
        new_rows[i].resize(matrix.size() - 1);
    }

    size_t i_new = 0, j_new = 0;
    for(size_t i = 0; i < matrix.size(); ++i) {
        if(i != e_row) {
            j_new = 0;
            for(size_t j = 0; j < matrix.size(); ++j) {
                if(j != e_column) {
                    new_rows[i_new][j_new] = matrix(i, j);
                    ++j_new;
                }
            }
            ++i_new;
        }
    }

    return Square_Matrix<long long int>(new_rows);
}

//naive determinant
//write with templates too lazily
long long int bad_determinant(const Square_Matrix<long long int>& matrix) {
    if(matrix.size() == 0) {
        return 0;
    }
    if(matrix.size() == 1) {
        return matrix(0, 0);
    }

    long long int det = 0;
    for(size_t i = 0; i < matrix.size(); ++i) {
        Square_Matrix<long long int> add_matrix = make_add_matrix(matrix, 0, i);
        long long int minor = bad_determinant(add_matrix);

        if(i % 2 == 1) {
            minor *= -1;
        }

        det += (minor * matrix(0, i));
    }

    return det;
}

Square_Matrix<long long int> matrix_generator(size_t matrix_size) {
    std::vector<std::vector<long long int>> rows;
    rows.resize(matrix_size);
    for(size_t i = 0; i < (matrix_size); ++i) {
        rows[i].resize(matrix_size);
    }

    for(size_t i = 0; i < matrix_size; ++i) {
        for(size_t j = 0; j < matrix_size; ++j) {
            rows[i][j] = (rand() % 10);
        }
    }

    return Square_Matrix<long long int>(rows);
}
}

void long_test(int n_tests = 1000) {
    srand(0);
    size_t sum_bad_time = 0, sum_norm_time = 0;

    for(int i = 0; i < n_tests; ++i) {
        Square_Matrix<long long int>m = matrix_generator(7);

        long start_time = clock();
        long long int bad_det = bad_determinant(m);
        long bad_time = clock() - start_time;

        start_time = clock();
        long long int det = determinant(m);
        long norm_time = clock() - start_time;

        assert(det == bad_det);
        std::cout << "Bad det = " << bad_det << std::endl;
        std::cout << "Det = " << det << std::endl;
        std::cout << "Bad time = " << bad_time << std::endl;
        std::cout << "Normal time = " << norm_time << std::endl;

        sum_bad_time += static_cast<size_t>(bad_time);
        sum_norm_time += static_cast<size_t>(norm_time);
    }

    std::cout << "Sum bad time = " << sum_bad_time << std::endl;
    std::cout << "Sum normal time = " << sum_norm_time << std::endl;
}
