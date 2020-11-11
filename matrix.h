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
    Square_Matrix(const T* inp_data, size_t inp_size);

    Square_Matrix(const std::vector<std::vector<T>>& input_rows);

    ~Square_Matrix();

    Square_Matrix(const Square_Matrix& matr): Square_Matrix(matr.data_, matr.size_) {}

    Square_Matrix& operator= (const Square_Matrix& matr)&;

    template<typename U> Square_Matrix(const Square_Matrix<U>& matr);

    void transpose() const;

    void add_row_to_row(size_t src_num, size_t dst_num) const;

    void sub_row_to_row(size_t src_num, size_t dst_num) const;

    T operator()(size_t row_i, size_t column_i) const;

    size_t size() const { return size_; }

    bool operator==(const Square_Matrix<T>& rhs) const;

    Square_Matrix& operator+=(const Square_Matrix& rhs)&;

    Square_Matrix& operator-=(const Square_Matrix<T>& rhs)&;

    Square_Matrix& operator*=(const Square_Matrix<T>& rhs)&;

    void print() const;

    T determinant() const;
};

//----------------------------methods for Square_Matrix--------------------------------
template<typename T>
Square_Matrix<T>::Square_Matrix(const T* inp_data, size_t inp_size):
    size_(inp_size)
{
    if(size_ > 0) data_ = new T[size_ * size_];

    for(size_t i = 0; i < size_ * size_; ++i) {
        data_[i] = inp_data[i];
    }
}

template<typename T>
Square_Matrix<T>::Square_Matrix(const std::vector<std::vector<T>>& input_rows):
    size_(input_rows.size())
{
    for(const std::vector<T> &row : input_rows) {
        assert((row.size() == size_) && ("Invalid matrix size"));
    }

    if(size_ > 0) data_ = new T[size_ * size_];

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i * size_ + j] = input_rows[i][j];
        }
    }
}

template<typename T>
Square_Matrix<T>::~Square_Matrix() {
    if(size_ > 0) delete [] data_;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator= (const Square_Matrix<T>& matr)& {
    if(this != &matr) {
        if((size_ > 0)) delete [] data_;

        size_ = matr.size();

        if((size_ > 0)) data_ = new T[size_ * size_];

        for(size_t i = 0; i < size_; ++i) {
            for(size_t j = 0; j < size_; ++j) {
                data_[i * size_ + j] = matr(i, j);
            }
        }
    }

    return *this;
}

template<typename T>
template<typename U>
Square_Matrix<T>::Square_Matrix(const Square_Matrix<U>& matr):
    size_(matr.size())
{
    if(size_ > 0) data_ = new T[size_ * size_];
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i * size_ + j] = static_cast<T>(matr(i, j));
        }
    }
}

template<typename T>
void Square_Matrix<T>::transpose() const {
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < i; ++j) {
            std::swap(data_[i * size_ + j], data_[j * size_ + i]);
        }
    }
}

template<typename T>
void Square_Matrix<T>::add_row_to_row(size_t src_num, size_t dst_num) const {
    assert((src_num < size_) && "invalid row number");
    assert((dst_num < size_) && "invalid row number");

    for(size_t i = 0; i < size_; ++i) {
        data_[dst_num * size_ + i] += data_[src_num * size_ + i];
    }
}

template<typename T>
void Square_Matrix<T>::sub_row_to_row(size_t src_num, size_t dst_num) const {
    assert((src_num < size_) && "invalid row number");
    assert((dst_num < size_) && "invalid row number");

    for(size_t i = 0; i < size_; ++i) {
        data_[dst_num * size_ + i] -= data_[src_num * size_ + i];
    }
}

template<typename T>
T Square_Matrix<T>::operator()(size_t row_i, size_t column_i) const {
    assert((row_i < size_) && "invalid row");
    assert((column_i < size_) && "invalid column");
    return data_[row_i * size_ + column_i];
}

template<typename T>
bool Square_Matrix<T>::operator==(const Square_Matrix<T>& rhs) const {

    if(size_ != rhs.size()) return false;

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            if(!is_match(data_[i * size_ + j], rhs(i, j))) return false;
        }
    }

    return true;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator+=(const Square_Matrix<T>& rhs)& {
    assert((size_ == rhs.size()) && "different matrix sizes");
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i * size_ + j] += rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator-=(const Square_Matrix<T>& rhs)& {
    assert((size_ == rhs.size()) && "different matrix sizes");
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i * size_ + j] -= rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator*=(const Square_Matrix<T>& rhs)& { //naive algorithm (almost)
    assert((size_ == rhs.size()) && "different matrix sizes");
    if(size_ == 0) return *this;

    T *new_data = new T[size_ * size_];

    Square_Matrix new_rhs(rhs);
    new_rhs.transpose(); //for cache friendly

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            new_data[i * size_ + j] = 0;
            for(size_t k = 0; k < size_; ++k) {
                new_data[i * size_ + j] += (data_[i * size_ + k] * new_rhs(j, k));
            }
        }
    }

    delete [] data_;
    data_ = new_data;

    return *this;
}

template<typename T>
void Square_Matrix<T>::print() const {
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            std::cout << data_[i * size_ + j] << " ";
        }
        std::cout << std::endl;
    }
}

template<typename T>
Square_Matrix<T> operator+(const Square_Matrix<T> lhs, const Square_Matrix<T> rhs) {
    Square_Matrix<T> tmp(lhs);
    tmp += rhs;
    return tmp;
}

template<typename T>
Square_Matrix<T> operator-(const Square_Matrix<T> lhs, const Square_Matrix<T> rhs) {
    Square_Matrix<T> tmp(lhs);
    tmp -= rhs;
    return tmp;
}

template<typename T>
Square_Matrix<T> operator*(const Square_Matrix<T> lhs, const Square_Matrix<T> rhs) {
    Square_Matrix<T> tmp(lhs);
    tmp *= rhs;
    return tmp;
}

template<>
double Square_Matrix<double>::determinant() const;

template<>
int Square_Matrix<int>::determinant() const;

template<>
long long int Square_Matrix<long long int>::determinant() const;

}
