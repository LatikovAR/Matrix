#pragma once

#include <string>
#include <fstream>
#include <ctime>
#include <cstdlib>

#include "matrix.h"

namespace matrix {

class Matrix_Chain_Generator final {
    const size_t MAX_MATRIX_SIZE = 50;
    const int MAX_MATRIX_ELEM = 10;

    void write_matr_to_file(std::ofstream& out, Matrix<int> m);
    Matrix<int> gen_matr(size_t column_size, size_t row_size);
public:
    Matrix_Chain_Generator(std::string filename, size_t chain_size);
};

Matrix_Chain_Generator::Matrix_Chain_Generator(std::string filename, size_t chain_size) {
    std::srand(std::time(0));

    std::ofstream out;
    out.open(filename);

    out << chain_size << std::endl;

    size_t column_size, row_size;
    column_size = std::rand() % MAX_MATRIX_SIZE;

    for(size_t i = 0; i < chain_size; ++i) {
        row_size = std::rand() % MAX_MATRIX_SIZE;

        Matrix<int> m = gen_matr(column_size, row_size);
        write_matr_to_file(out, std::move(m));

        column_size = row_size;
    }

    out.close();
}

void Matrix_Chain_Generator::write_matr_to_file(std::ofstream& out, Matrix<int> m) {
    out << m.column_size() << " " << m.row_size() << std::endl;
    for(size_t i = 0; i < m.column_size(); ++i) {
        for(size_t j = 0; j < m.row_size(); ++j) {
            out << m(i, j) << " ";
        }
        out << std::endl;
    }
}

Matrix<int> Matrix_Chain_Generator::gen_matr(size_t column_size, size_t row_size) {
    Matrix<int> m(column_size, row_size);

    for(size_t i = 0; i < column_size; ++i) {
        for(size_t j = 0; j < row_size; ++j) {
            m(i, j) = std::rand() % MAX_MATRIX_ELEM;
        }
    }

    return m;
}

} //namespace matrix
