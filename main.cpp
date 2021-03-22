#include <iostream>
#include <cstdlib>
#include <chrono>

#include "matrix/matrix.h"
#include "matrix/matrix_mult_chain.h"
#include "matrix/matrix_chain_test.h"
#include "matrix/matrix_chain_generator.h"

using namespace matrix;

Matrix<int> input_matrix() {
    size_t column_size, row_size;
    std::cin >> column_size >> row_size;
    Matrix<int> m(column_size, row_size);

    for(size_t i = 0; i < column_size; ++i) {
        for(size_t j = 0; j < row_size; ++j) {
            std::cin >> m(i, j);
        }
    }

    return m;
}

int main() {
    //Matrix_Chain_Generator("test1.txt", 5);

    Matrix_Chain<int> chain;
    size_t chain_size;
    std::cin >> chain_size;

    for(size_t i = 0; i < chain_size; ++i) {
        chain.insert_matrix(input_matrix(), i);
    }

    //optimal way
    auto start_time = std::chrono::high_resolution_clock::now();

    auto optimal_trace = chain.compute_optimal_trace();
    Matrix<int> m1 = chain.mult_chain_optimal(optimal_trace);

    auto end_time = std::chrono::high_resolution_clock::now();

    double optimal_time =
            std::chrono::duration<double,
            std::chrono::milliseconds::period>
            (end_time - start_time).count();

    //naive way
    start_time = std::chrono::high_resolution_clock::now();

    Matrix<int> m2 = chain.mult_chain_naive();

    end_time = std::chrono::high_resolution_clock::now();

    double naive_time =
            std::chrono::duration<double,
            std::chrono::milliseconds::period>
            (end_time - start_time).count();

    assert(m1 == m2);
    std::cout << "optimal time = " << optimal_time << std::endl;
    std::cout << "naive time = " << naive_time << std::endl;

    Matrix_Chain<int>::print_optimal_trace_with_brackets(optimal_trace);
    return 0;
}
