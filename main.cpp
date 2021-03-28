#include <iostream>
#include <cstdlib>
#include <chrono>

#include "matrix/matrix.h"
#include "matrix/matrix_mult_chain.h"
#include "matrix/matrix_chain_test.h"
#include "matrix/matrix_chain_generator.h"

using namespace matrix;

/*
int main() {
    matrix_test::add_test();
    matrix_test::trace_test();
    return 0;
}*/

void gen_numbers_for_matrix(Matrix<int>& m) {
    const int MAX_NUM = 10;
    for(size_t i = 0; i < m.column_size(); ++i) {
        for(size_t j = 0; j < m.row_size(); ++j) {
            m(i, j) = rand() % MAX_NUM;
        }
    }
}

int main() {
    Matrix_Chain<int> chain;
    srand(time(0));

    while(!std::cin.eof()) {
        size_t column_size, row_size;
        std::cin >> column_size >> row_size;
        if(std::cin.eof()) break;

        Matrix<int> m(column_size, row_size);
        gen_numbers_for_matrix(m);

        try {
            chain.insert_matrix(std::move(m), chain.matrix_num());
        }  catch (std::invalid_argument& e) {
            std::cerr << e.what();
            exit(-1);
        }
    }

    //optimal way
    auto start_time = std::chrono::high_resolution_clock::now();

    auto optimal_trace = chain.compute_optimal_trace();
    //start_time = std::chrono::high_resolution_clock::now();
    Matrix<int> m1 = chain.mult_chain_optimal(optimal_trace.op_order());

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
    std::cout << "chain size = " << chain.matrix_num() << std::endl;
    std::cout << "optimal time = " << optimal_time << std::endl;
    std::cout << "naive time = " << naive_time << std::endl;

    std::cout << "optimal * number = " << optimal_trace.weight << std::endl;
    std::cout << "naive * number = " << chain.naive_weight() << std::endl;

    std::cout << "optimal multiplication order: ";
    Matrix_Chain<int>::print_optimal_trace_with_numbers(optimal_trace.op_order());
    return 0;
}
