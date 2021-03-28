#include <cassert>
#include <stdexcept>
#include <iostream>

#include "matrix.h"
#include "matrix_mult_chain.h"
#include "matrix_chain_test.h"

using namespace matrix;

namespace matrix_test {

void add_test() {
    Matrix_Chain<int> matrix_chain;

    Matrix<int> m1(2, 2);
    m1(0, 0) = 1;
    m1(0, 1) = 1;
    m1(1, 0) = 2;
    m1(1, 1) = 2;

    matrix_chain.insert_matrix(m1, 0);

    assert(m1 == matrix_chain.ret_matrix(0));

    Matrix<int> m2(2, 3);
    m2(0, 0) = 3;
    m2(0, 1) = 3;
    m2(0, 2) = 3;
    m2(1, 0) = 2;
    m2(1, 1) = 2;
    m2(1, 1) = 2;

    matrix_chain.insert_matrix(m2, 1);
    assert(m2 == matrix_chain.ret_matrix(1));
    assert(m1 == matrix_chain.ret_matrix(0));

    try {
        matrix_chain.insert_matrix(m2, 1);
        assert(0);
    }  catch (std::invalid_argument& e) {
        assert(e.what() == std::string("Matrix_chain: invalid inserting matrix size"));
    }

    try {
        matrix_chain.insert_matrix(m2, 0);
        assert(0);
    }  catch (std::invalid_argument& e) {
        assert(e.what() == std::string("Matrix_chain: invalid inserting matrix size"));
    }

    std::cout << "Add test complete\n";
}

void trace_test() {
    Matrix<int> m1(2, 2);
    Matrix<int> m2(2, 8);
    Matrix<int> m3(8, 3);
    Matrix<int> m4(3, 4);

    Matrix_Chain<int>  matrix_chain;
    matrix_chain.insert_matrix(m4, 0);
    matrix_chain.insert_matrix(m3, 0);
    matrix_chain.insert_matrix(m2, 0);
    matrix_chain.insert_matrix(m1, 0);

    std::vector<size_t> order = matrix_chain.compute_optimal_trace().op_order();
    Matrix_Chain<int>::print_optimal_trace_with_numbers(order);
    Matrix_Chain<int>::print_optimal_trace_with_brackets(order);
}

} //namespace matrix
