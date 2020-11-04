#pragma once

#include <cstdlib>
#include <vector>
#include <cassert>
#include <utility>

namespace matrix {

template <typename T> class Matrix {
protected:
    std::vector<std::vector<T>> rows;
    size_t column_size;
    size_t row_size;
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

    void transpose() {
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
        rows = new_rows;

        std::swap(row_size, column_size);
    }

    void add_row_to_row(size_t src_num, size_t dst_num) {
        assert((src_num < column_size) && "invalid row number");
        assert((dst_num < column_size) && "invalid row number");

        for(size_t i = 0; i < row_size; ++i) {
            rows[dst_num][i] += rows[src_num][i];
        }
    }

    T operator()(size_t row_i, size_t column_i) const {
        assert((row_i < column_size) && "invalid row");
        assert((column_i < row_size) && "invalid column");
        return rows[row_i][column_i];
    }

    size_t ret_column_size() const { return column_size; }
    size_t ret_row_size() const { return row_size; }

    bool operator==(Matrix m) {
        return (row_size == m.row_size) && (column_size == m.column_size) && (rows == m.rows);
    }
};

template <typename T> class Square_Matrix final: public Matrix<T> {
public:
    Square_Matrix(const std::vector<std::vector<T>>& input_rows): Matrix<T>(input_rows) {
        assert((Matrix<T>::column_size == Matrix<T>::row_size) && "It isn't square matrix");
    }

    ~Square_Matrix() {}
};

double determinant(const Square_Matrix<double>& matrix);

}
