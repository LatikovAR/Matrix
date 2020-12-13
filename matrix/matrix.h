#pragma once

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <utility>
#include <cmath>
#include <typeinfo>
#include <stdexcept>

#include "../mem_storage/mem_storage.h"

namespace matrix {

template <typename T> class Abstract_Matrix {
protected:
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
    virtual ~Abstract_Matrix() {}
    virtual const T& operator()(size_t row_i, size_t column_i) const = 0;
    virtual bool operator==(const Abstract_Matrix<T>& rhs) const = 0;
    virtual void print() const = 0;
};

template<typename T>
class Symmetric_Matrix;

template <typename T>
class Square_Matrix final :
        public Abstract_Matrix<T>,
        protected mem_storage::Matrix_Storage<T>
{
private:
    const size_t& size_ = mem_storage::Matrix_Storage<T>::column_size_;
    using mem_storage::Matrix_Storage<T>::data_;
    using mem_storage::Matrix_Storage<T>::used_;
    using mem_storage::Matrix_Storage<T>::swap;

    using Abstract_Matrix<T>::is_match;

public:
    Square_Matrix(const T *const *data, size_t size);

    Square_Matrix(const std::vector<std::vector<T>>& input_rows);

    Square_Matrix(const Square_Matrix& rhs);

    Square_Matrix(const Symmetric_Matrix<T>& rhs);

    Square_Matrix(size_t size = 0);

    template<typename U> Square_Matrix(const Square_Matrix<U>& rhs);

    Square_Matrix& operator= (const Square_Matrix& rhs)&;

    void transpose() const;

    void add_row_to_row(size_t src_num, size_t dst_num) const;

    void sub_row_to_row(size_t src_num, size_t dst_num, T koef = static_cast<T>(1)) const;

    void swap_rows(size_t row1_num, size_t row2_num) const;

    void mult_row_on_number(size_t row_num, T number) const;

    const T& operator()(size_t row_i, size_t column_i) const override;

    T& operator()(size_t row_i, size_t column_i);

    size_t size() const { return size_; }

    bool operator==(const Abstract_Matrix<T>& inp_rhs) const override;

    Square_Matrix& operator+=(const Square_Matrix& rhs)&;

    Square_Matrix& operator-=(const Square_Matrix<T>& rhs)&;

    Square_Matrix& operator*=(const Square_Matrix<T>& rhs)&;

    void print() const override;

    T determinant() const;
};


template <typename T>
class Symmetric_Matrix final :
        public Abstract_Matrix<T>,
        protected mem_storage::Symmetric_Matrix_Storage<T>
{
private:
    const size_t& size_ = mem_storage::Symmetric_Matrix_Storage<T>::size_;
    using mem_storage::Symmetric_Matrix_Storage<T>::data_;
    using mem_storage::Symmetric_Matrix_Storage<T>::used_;
    using mem_storage::Symmetric_Matrix_Storage<T>::swap;

    using Abstract_Matrix<T>::is_match;

public:
    Symmetric_Matrix(size_t size = 0);

    Symmetric_Matrix(const Symmetric_Matrix& rhs);

    Symmetric_Matrix& operator= (const Symmetric_Matrix& rhs)&;

    template<typename U> Symmetric_Matrix(const Symmetric_Matrix<U>& rhs);

    const T& operator()(size_t row_i, size_t column_i) const override;

    T& operator()(size_t row_i, size_t column_i);

    size_t size() const { return size_; }

    bool operator==(const Abstract_Matrix<T>& inp_rhs) const override;

    Symmetric_Matrix& operator+=(const Symmetric_Matrix& rhs)&;

    Symmetric_Matrix& operator-=(const Symmetric_Matrix<T>& rhs)&;

    void print() const override;

    T determinant() const;
};


template <typename T>
class Matrix final :
        public Abstract_Matrix<T>,
        protected mem_storage::Matrix_Storage<T>
{
private:
    using mem_storage::Matrix_Storage<T>::column_size_;
    using mem_storage::Matrix_Storage<T>::row_size_;
    using mem_storage::Matrix_Storage<T>::data_;
    using mem_storage::Matrix_Storage<T>::used_;
    using mem_storage::Matrix_Storage<T>::swap;

    using Abstract_Matrix<T>::is_match;
public:
    Matrix() {}

    //matrix elems should be contained in inp_data row by row
    Matrix(const T *const *data, size_t column_size, size_t row_size);

    Matrix(const std::vector<std::vector<T>>& input_rows);

    Matrix(const Matrix& rhs);

    Matrix(size_t column_size, size_t row_size);

    Matrix& operator= (const Matrix& rhs)&;

    Matrix(const Square_Matrix<T>& rhs);

    Matrix(const Symmetric_Matrix<T>& rhs);

    template<typename U> Matrix(const Matrix<U>& rhs);

    void transpose();

    void add_row_to_row(size_t src_num, size_t dst_num) const;

    void sub_row_to_row(size_t src_num, size_t dst_num, T koef = static_cast<T>(1)) const;

    void swap_rows(size_t row1_num, size_t row2_num) const;

    void mult_row_on_number(size_t row_num, T number) const;

    const T& operator()(size_t row_i, size_t column_i) const override;

    T& operator()(size_t row_i, size_t column_i);

    size_t row_size() const { return row_size_; }
    size_t column_size() const { return column_size_; }

    bool operator==(const Abstract_Matrix<T>& inp_rhs) const override;

    Matrix& operator+=(const Matrix<T>& rhs)&;

    Matrix& operator-=(const Matrix<T>& rhs)&;

    Matrix& operator*=(const Matrix<T>& rhs)&;

    void print() const override;
};


class Linear_Equations_System final {
private:
    //linear system
    const Matrix<double> matrix_;
    const std::vector<double> column_;

    Matrix<double> matrix_buf_;
    std::vector<double> column_buf_;
    std::vector<bool> trace_; //for going back the same way as straight pass
    size_t row_passed_;
    size_t column_passed_;
    bool has_solution_ = true;
    std::vector<double> answer_;

    static constexpr double DOUBLE_GAP = 1e-12;

    static bool is_match(double a, double b) { //instead of ==
        return (fabs(a - b) < DOUBLE_GAP);
    }

    void straight_pass();
    void reverse_pass();
public:
    Linear_Equations_System(const Matrix<double>& matrix, const std::vector<double>& column):
        matrix_(matrix), column_(column),
        matrix_buf_(matrix), column_buf_(column)
    {
        assert((matrix.column_size() == column.size()) && "invalid data sizes");
    }

    //returns any of the solutions if it isn't one
    //if no solution returning result undefined and bool = false
    //class object condition undefined after this method
    std::pair<std::vector<double>, bool> solve();


};


//----------------------------methods for Square_Matrix--------------------------------

template<typename T>
Square_Matrix<T>::Square_Matrix(const T *const *data, size_t size):
    mem_storage::Matrix_Storage<T>(size, size)
{
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            new (&(data_[i][j])) T(data[i][j]);
            ++used_;
        }
    }
}

template<typename T>
Square_Matrix<T>::Square_Matrix(const std::vector<std::vector<T>>& input_rows):
    mem_storage::Matrix_Storage<T>(input_rows.size(), input_rows.size())
{
    for(const std::vector<T> &row : input_rows) {
        if(row.size() != size_)
            throw std::invalid_argument
                ("Constructing matrix from non-suitable std::vector<std::vector<T>>");
    }

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            new (&(data_[i][j])) T(input_rows[i][j]);
            ++used_;
        }
    }
}

template<typename T>
Square_Matrix<T>::Square_Matrix(const Square_Matrix& rhs):
    mem_storage::Matrix_Storage<T>(rhs.size(), rhs.size())
{
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            new (&(data_[i][j])) T(rhs(i, j));
            ++used_;
        }
    }
}

template<typename T>
Square_Matrix<T>::Square_Matrix(const Symmetric_Matrix<T>& rhs):
    mem_storage::Matrix_Storage<T>(rhs.size(), rhs.size())
{
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            new (&(data_[i][j])) T(rhs(i, j));
            ++used_;
        }
    }
}

template<typename T>
Square_Matrix<T>::Square_Matrix(size_t size):
    mem_storage::Matrix_Storage<T>(size, size)
{
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
           new (&(data_[i][j])) T();
            ++used_;
        }
    }
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator= (const Square_Matrix<T>& rhs)& {
    Square_Matrix<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T>
template<typename U>
Square_Matrix<T>::Square_Matrix(const Square_Matrix<U>& rhs):
    mem_storage::Matrix_Storage<T>(rhs.size(), rhs.size())
{
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            new (&(data_[i][j])) T(static_cast<T>(rhs(i, j)));
            ++used_;
        }
    }
}

template<typename T>
void Square_Matrix<T>::transpose() const {
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < i; ++j) {
            std::swap(data_[i][j], data_[j][i]);
        }
    }
}

template<typename T>
void Square_Matrix<T>::add_row_to_row(size_t src_num, size_t dst_num) const {
    if((src_num >= size_) || (dst_num >= size_))
        throw std::out_of_range("invalid matrix row number");

    for(size_t i = 0; i < size_; ++i) {
        data_[dst_num][i] += data_[src_num][i];
    }
}

template<typename T>
void Square_Matrix<T>::sub_row_to_row(size_t src_num, size_t dst_num, T koef) const {
    if((src_num >= size_) || (dst_num >= size_))
        throw std::out_of_range("invalid matrix row number");

    for(size_t i = 0; i < size_; ++i) {
        data_[dst_num][i] -= (data_[src_num][i] * koef);
    }
}

template<typename T>
void Square_Matrix<T>::swap_rows(size_t row1_num, size_t row2_num) const {
    if((row1_num >= size_) || (row2_num >= size_))
        throw std::out_of_range("invalid matrix row number");
    std::swap(data_[row1_num], data_[row2_num]);
}

template<typename T>
void Square_Matrix<T>::mult_row_on_number(size_t row_num, T number) const {
    if(row_num >= size_)
        throw std::out_of_range("invalid matrix row number");
    for(size_t i = 0; i < size_; ++i) {
        data_[row_num][i] *= number;
    }
}

template<typename T>
const T& Square_Matrix<T>::operator()(size_t row_num, size_t column_num) const {
    if((row_num >= size_) || (column_num >= size_))
        throw std::out_of_range("invalid matrix elem");
    return data_[row_num][column_num];
}

template<typename T>
T& Square_Matrix<T>::operator()(size_t row_num, size_t column_num) {
    if((row_num >= size_) || (column_num >= size_))
        throw std::out_of_range("invalid matrix elem");
    return data_[row_num][column_num];
}

template<typename T>
bool Square_Matrix<T>::operator==(const Abstract_Matrix<T>& inp_rhs) const {
    if(typeid(*this) != typeid(inp_rhs)) return false;

    const Square_Matrix<T>& rhs = static_cast<const Square_Matrix&>(inp_rhs);

    if(size_ != rhs.size()) return false;

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            if(!is_match(data_[i][j], rhs(i, j))) return false;
        }
    }

    return true;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator+=(const Square_Matrix<T>& rhs)& {
    if(size_ != rhs.size_)
        throw std::invalid_argument("invalid matrices' sizes for +");
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i][j] += rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator-=(const Square_Matrix<T>& rhs)& {
    if(size_ != rhs.size_)
        throw std::invalid_argument("invalid matrices' sizes for -");
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i][j] -= rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator*=(const Square_Matrix<T>& rhs)& { //naive algorithm (almost)
    if(size_ != rhs.size_)
        throw std::invalid_argument("invalid matrices' sizes for *");
    if(size_ == 0) return *this;

    Square_Matrix<T> res_tmp(size_);
    Square_Matrix rhs_copy(rhs);
    rhs_copy.transpose(); //for cache friendly

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            res_tmp(i, j) = 0;
            for(size_t k = 0; k < size_; ++k) {
                res_tmp(i, j) += (data_[i][k] * rhs_copy(j, k));
            }
        }
    }

    swap(res_tmp);
    return *this;
}

template<typename T>
void Square_Matrix<T>::print() const {
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            std::cout << data_[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

template<typename T>
Square_Matrix<T> operator+(const Square_Matrix<T>& lhs, const Square_Matrix<T>& rhs) {
    Square_Matrix<T> tmp(lhs);
    tmp += rhs;
    return tmp;
}

template<typename T>
Square_Matrix<T> operator-(const Square_Matrix<T>& lhs, const Square_Matrix<T>& rhs) {
    Square_Matrix<T> tmp(lhs);
    tmp -= rhs;
    return tmp;
}

template<typename T>
Square_Matrix<T> operator*(const Square_Matrix<T>& lhs, const Square_Matrix<T>& rhs) {
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




//----------------------------methods for Symmetric_Matrix--------------------------------

template<typename T>
Symmetric_Matrix<T>::Symmetric_Matrix(const Symmetric_Matrix& rhs):
    mem_storage::Symmetric_Matrix_Storage<T>(rhs.size())
{
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < (i + 1); ++j) {
            new(&(data_[i][j])) T(rhs(i, j));
            ++used_;
        }
    }
}

template<typename T>
Symmetric_Matrix<T>::Symmetric_Matrix(size_t size):
    mem_storage::Symmetric_Matrix_Storage<T>(size)
{
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < (i + 1); ++j) {
            new(&(data_[i][j])) T();
            ++used_;
        }
    }
}

template<typename T>
Symmetric_Matrix<T>& Symmetric_Matrix<T>::operator= (const Symmetric_Matrix<T>& rhs)& {
    Symmetric_Matrix<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T>
template<typename U>
Symmetric_Matrix<T>::Symmetric_Matrix(const Symmetric_Matrix<U>& rhs):
    mem_storage::Symmetric_Matrix_Storage<T>(rhs.size())
{
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < (i + 1); ++j) {
            new (&(data_[i][j])) T(static_cast<T>(rhs(i, j)));
            ++used_;
        }
    }
}

template<typename T>
const T& Symmetric_Matrix<T>::operator()(size_t row_num, size_t column_num) const {
    if((row_num >= size_) || (column_num >= size_))
        throw std::out_of_range("invalid matrix elem");

    if(row_num >= column_num) return data_[row_num][column_num];
    return data_[column_num][row_num];
}

template<typename T>
T& Symmetric_Matrix<T>::operator()(size_t row_num, size_t column_num) {
    if((row_num >= size_) || (column_num >= size_))
        throw std::out_of_range("invalid matrix elem");

    if(row_num >= column_num) return data_[row_num][column_num];
    return data_[column_num][row_num];
}

template<typename T>
bool Symmetric_Matrix<T>::operator==(const Abstract_Matrix<T>& inp_rhs) const {
    if(typeid(*this) != typeid(inp_rhs)) return false;

    const Square_Matrix<T>& rhs = static_cast<const Symmetric_Matrix&>(inp_rhs);

    if(size_ != rhs.size()) return false;

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < (i + 1); ++j) {
            if(!is_match(data_[i][j], rhs(i, j))) return false;
        }
    }

    return true;
}

template<typename T>
Symmetric_Matrix<T>& Symmetric_Matrix<T>::operator+=(const Symmetric_Matrix<T>& rhs)& {
    if(size_ != rhs.size_)
        throw std::invalid_argument("invalid matrices' sizes for +");

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < (i + 1); ++j) {
            data_[i][j] += rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Symmetric_Matrix<T>& Symmetric_Matrix<T>::operator-=(const Symmetric_Matrix<T>& rhs)& {
    if(size_ != rhs.size_)
        throw std::invalid_argument("invalid matrices' sizes for -");

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < (i + 1); ++j) {
            data_[i][j] -= rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
void Symmetric_Matrix<T>::print() const {
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            std::cout << operator()(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

template<typename T>
Symmetric_Matrix<T> operator+(const Symmetric_Matrix<T>& lhs,
                              const Symmetric_Matrix<T>& rhs)
{
    Symmetric_Matrix<T> tmp(lhs);
    tmp += rhs;
    return tmp;
}

template<typename T>
Symmetric_Matrix<T> operator-(const Symmetric_Matrix<T>& lhs,
                              const Symmetric_Matrix<T>& rhs)
{
    Symmetric_Matrix<T> tmp(lhs);
    tmp -= rhs;
    return tmp;
}

template<>
double Symmetric_Matrix<double>::determinant() const;

template<>
int Symmetric_Matrix<int>::determinant() const;

template<>
long long int Symmetric_Matrix<long long int>::determinant() const;






//----------------------------methods for Matrix----------------------------------------

template<typename T>
Matrix<T>::Matrix(const T *const *data, size_t column_size, size_t row_size):
    mem_storage::Matrix_Storage<T>(column_size, row_size)
{
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            new (&(data_[i][j])) T(data[i][j]);
            ++used_;
        }
    }
}

template<typename T>
Matrix<T>::Matrix(const std::vector<std::vector<T>>& input_rows):
    mem_storage::Matrix_Storage<T>(
        input_rows.size(),
        (input_rows.size() > 0) ? input_rows[0].size() : 0)
{
    for(const std::vector<T> &row : input_rows) {
        if(row.size() != row_size_)
            throw std::invalid_argument("Constructing matrix from non-suitable std::vector<std::vector<T>>");
    }

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            new (&(data_[i][j])) T(input_rows[i][j]);
            ++used_;
        }
    }
}

template<typename T>
Matrix<T>::Matrix(const Matrix& rhs):
    mem_storage::Matrix_Storage<T>(rhs.column_size(), rhs.row_size())
{
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            new (&(data_[i][j])) T(rhs(i, j));
            ++used_;
        }
    }
}

template<typename T>
Matrix<T>::Matrix(size_t column_size, size_t row_size):
    mem_storage::Matrix_Storage<T>(column_size, row_size)
{
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            new (&(data_[i][j])) T();
            ++used_;
        }
    }
}


template<typename T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& rhs)& {
    Matrix<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T>
template<typename U>
Matrix<T>::Matrix(const Matrix<U>& rhs):
    mem_storage::Matrix_Storage<T>(rhs.column_size(), rhs.row_size())
{
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            new (&(data_[i][j])) T(static_cast<T>(rhs(i, j)));
            ++used_;
        }
    }
}

template<typename T>
Matrix<T>::Matrix(const Square_Matrix<T>& rhs):
    mem_storage::Matrix_Storage<T>(rhs.size(), rhs.size())
{
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            new (&(data_[i][j])) T(rhs(i, j));
            ++used_;
        }
    }
}

template<typename T>
Matrix<T>::Matrix(const Symmetric_Matrix<T>& rhs):
    mem_storage::Matrix_Storage<T>(rhs.size(), rhs.size())
{
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            new (&(data_[i][j])) T(rhs(i, j));
            ++used_;
        }
    }
}

template<typename T>
void Matrix<T>::transpose() {
    Matrix<T> tmp(row_size_, column_size_);

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            tmp(j, i) = data_[i][j];
        }
    }

    swap(tmp);
}

template<typename T>
void Matrix<T>::add_row_to_row(size_t src_num, size_t dst_num) const {
    if((src_num >= column_size_) || (dst_num >= column_size_))
        throw std::out_of_range("invalid matrix row number");

    for(size_t i = 0; i < row_size_; ++i) {
        data_[dst_num][i] += data_[src_num][i];
    }
}

template<typename T>
void Matrix<T>::sub_row_to_row(size_t src_num, size_t dst_num, T koef) const {
    if((src_num >= column_size_) || (dst_num >= column_size_))
        throw std::out_of_range("invalid matrix row number");

    for(size_t i = 0; i < row_size_; ++i) {
        data_[dst_num][i] -= (data_[src_num][i] * koef);
    }
}

template<typename T>
void Matrix<T>::swap_rows(size_t row1_num, size_t row2_num) const {
    if((row1_num >= column_size_) || (row2_num >= column_size_))
        throw std::out_of_range("invalid matrix row number");
    std::swap(data_[row1_num], data_[row2_num]);
}

template<typename T>
void Matrix<T>::mult_row_on_number(size_t row_num, T number) const {
    if(row_num >= column_size_)
        throw std::out_of_range("invalid matrix row number");
    for(size_t i = 0; i < row_size_; ++i) {
        data_[row_num][i] *= number;
    }
}

template<typename T>
const T& Matrix<T>::operator()(size_t column_i, size_t row_i) const {
    if((column_i >= column_size_) || (row_i >= row_size_))
        throw std::out_of_range("invalid matrix elem");
    return data_[column_i][row_i];
}

template<typename T>
T& Matrix<T>::operator()(size_t column_i, size_t row_i) {
    if((column_i >= column_size_) || (row_i >= row_size_))
        throw std::out_of_range("invalid matrix elem");
    return data_[column_i][row_i];
}

template<typename T>
bool Matrix<T>::operator==(const Abstract_Matrix<T>& inp_rhs) const {
    if(typeid(*this) != typeid(inp_rhs)) return false;

    const Matrix<T>& rhs = static_cast<const Matrix&>(inp_rhs);

    if((row_size_ != rhs.row_size()) || (column_size_ != rhs.column_size())) return false;

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            if(!is_match(data_[i][j], rhs(i, j))) return false;
        }
    }

    return true;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs)& {
    if((row_size_ != rhs.row_size()) || (column_size_ != rhs.column_size()))
        throw std::invalid_argument("invalid matrices' sizes for +");
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            data_[i][j] += rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs)& {
    if((row_size_ != rhs.row_size()) || (column_size_ != rhs.column_size()))
        throw std::invalid_argument("invalid matrices' sizes for -");
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            data_[i][j] -= rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& rhs)& {
    if(row_size_ != rhs.column_size())
        throw std::invalid_argument("invalid matrices' sizes for *");

    if((row_size_ == 0) || (column_size_ == 0)) return *this;

    if(rhs.row_size() == 0) {
        Matrix<T> tmp(column_size_, 0);
        swap(tmp);
        return *this;
    }

    Matrix<T> res_tmp(column_size_, rhs.row_size());
    Matrix<T> rhs_copy(rhs);
    rhs_copy.transpose(); //for cache friendly

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < rhs_copy.column_size(); ++j) {
            res_tmp(i, j) = 0;
            for(size_t k = 0; k < row_size_; ++k) {
                res_tmp(i, j) += (data_[i][k] * rhs_copy(j, k));
            }
        }
    }

    swap(res_tmp);
    return *this;
}

template<typename T>
void Matrix<T>::print() const {
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            std::cout << data_[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

template<typename T>
Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    Matrix<T> tmp(lhs);
    tmp += rhs;
    return tmp;
}

template<typename T>
Matrix<T> operator-(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    Matrix<T> tmp(lhs);
    tmp -= rhs;
    return tmp;
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    Matrix<T> tmp(lhs);
    tmp *= rhs;
    return tmp;
}

//----------------------------------------Solving_linear_systems--------------------------------

template<typename T1, typename T2>
std::vector<T2> vector_conversion(const std::vector<T1>& inp_vector) {
    std::vector<T2> res_vector;
    res_vector.resize(inp_vector.size());
    for(size_t i = 0; i < inp_vector.size(); ++i) {
        res_vector[i] = static_cast<T2>(inp_vector[i]);
    }
    return res_vector;
}

template<typename T>
std::pair<std::vector<double>, bool> solve_linear_equations(const Square_Matrix<T>& matrix,
                                                            const std::vector<T>& column) {
    return Linear_Equations_System(Matrix<double>(Square_Matrix<T>(matrix)),
                                   vector_conversion<T, double>(column)).solve();
}


template<typename T>
std::pair<std::vector<double>, bool> solve_linear_equations(const Matrix<T>& matrix,
                                                            const std::vector<T>& column) {
    return Linear_Equations_System(Matrix<double>(matrix),
                                   vector_conversion<T, double>(column)).solve();
}
} //namespace matrix
