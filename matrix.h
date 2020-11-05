#pragma once

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <utility>
#include <cmath>
#include <new>

namespace matrix {

template <typename T> class Square_Matrix final {
private:
    const size_t size_;
    T *data_;

    static constexpr double DOUBLE_GAP = 1e-12;

    static bool is_match(double lhs, double rhs) { //instead of ==
        return (fabs(lhs - rhs) < DOUBLE_GAP);
    }

    static bool is_match(int lhs, int rhs) {
        return (lhs == rhs);
    }

    static bool is_match(long long int lhs, long long int rhs) {
        return (lhs == rhs);
    }
public:
    Square_Matrix(const std::vector<std::vector<T>>& input_rows):
        size_(input_rows.size()),
        data_(new T[size_ * size_])
    {
        for(const std::vector<T> &row : input_rows) {
            assert((row.size() == size_) && ("Invalid matrix size"));
        }

        for(size_t i = 0; i < size_; ++i) {
            for(size_t j = 0; j < size_; ++j) {
                data_[i * size_ + j] = input_rows[i][j];
            }
        }
    }

    ~Square_Matrix() {
        delete [] data_;
    }

    void transpose() const {
        for(size_t i = 0; i < size_; ++i) {
            for(size_t j = 0; j < i; ++j) {
                std::swap(data_[i * size_ + j], data_[j * size_ + i]);
            }
        }
    }

    void add_row_to_row(size_t src_num, size_t dst_num) const {
        assert((src_num < size_) && "invalid row number");
        assert((dst_num < size_) && "invalid row number");

        for(size_t i = 0; i < size_; ++i) {
            data_[dst_num * size_ + i] += data_[src_num * size_ + i];
        }
    }

    T operator()(size_t row_i, size_t column_i) const {
        assert((row_i < size_) && "invalid row");
        assert((column_i < size_) && "invalid column");
        return data_[row_i * size_ + column_i];
    }

    size_t size() const { return size_; }

    template<typename U> Square_Matrix(const Square_Matrix<U>& matr): size_(matr.size()) {
        data_ = new T[size_ * size_];
        for(size_t i = 0; i < size_; ++i) {
            for(size_t j = 0; j < size_; ++j) {
                data_[i * size_ + j] = static_cast<T>(matr(i, j));
            }
        }
    }

    bool operator==(const Square_Matrix<T>& matr) const {
        if(size_ != matr.size()) return false;

        for(size_t i = 0; i < size_; ++i) {
            for(size_t j = 0; j < size_; ++j) {
                if(!is_match(data_[i * size_ + j], matr(i, j))) return false;
            }
        }

        return true;
    }
};

double determinant(const Square_Matrix<double>& matrix);
float determinant(const Square_Matrix<float>& matrix);
int determinant(const Square_Matrix<int>& matrix);
long long int determinant(const Square_Matrix<long long int>& input_matrix);
}
