#pragma once

#include <cstdlib>
#include <vector>
#include <cassert>

namespace matrix {

template <typename T> class Matrix {
protected:
    const std::vector<std::vector<T>> rows;
    const size_t column_size;
    const size_t row_size;
public:
    Matrix(const std::vector<std::vector<T>>& input_rows):
        rows(input_rows),
        column_size(input_rows.size()),
        row_size((input_rows.size() > 0) ? input_rows[0].size() : 0)
    {
        for(const std::vector<T> &row : input_rows) {
            assert((row.size() == row_size) && ("Different row.size() in matrix isn't available"));
        }
    }

    virtual ~Matrix() {}

    Matrix transpose() const {
        std::vector<std::vector<T>> new_rows;
        std::vector<T> row_buffer;

        new_rows.resize(row_size);
        row_buffer.resize(column_size);

        for(size_t column_num = 0; column_num < row_size; ++column_num) {
            for(size_t row_num = 0; row_num < column_size; ++row_num) {
                row_buffer[row_num] = rows[row_num][column_num];
            }
            new_rows[column_num] = row_buffer;
        }

        Matrix new_matrix(new_rows);
        return new_matrix;
    }

    bool operator==(Matrix m) {
        return (row_size == m.row_size) && (column_size == m.column_size) && (rows == m.rows);
    }
};

template <typename T> class Square_Matrix final: public Matrix<T> {
public:
    Square_Matrix(const std::vector<std::vector<T>>& input_rows): Matrix<T>(input_rows) {
        assert(Matrix<T>::column_size == Matrix<T>::row_size);
    }

    ~Square_Matrix() {}


};

}
