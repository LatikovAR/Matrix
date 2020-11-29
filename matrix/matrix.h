#pragma once

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <utility>
#include <cmath>
#include <new>
#include <typeinfo>

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
    virtual T operator()(size_t row_i, size_t column_i) const = 0;
    virtual bool operator==(const Abstract_Matrix<T>& rhs) const = 0;
    virtual void print() const = 0;
};

template <typename T> class Matrix; //for conversion between different matrices

template <typename T>
class Square_Matrix final : public Abstract_Matrix<T> {
private:
    size_t size_ = 0;
    T **data_ = nullptr;

    using Abstract_Matrix<T>::is_match;

    friend class Matrix<T>;

    void new_data();

    static T** new_data_buf(size_t size);

    void delete_data();
public:
    Square_Matrix(const T *const *inp_data, size_t inp_size);

    Square_Matrix(const std::vector<std::vector<T>>& input_rows);

    ~Square_Matrix() override;

    Square_Matrix(const Square_Matrix& matr): Square_Matrix(matr.data_, matr.size_) {}

    Square_Matrix& operator= (const Square_Matrix& matr)&;

    template<typename U> Square_Matrix(const Square_Matrix<U>& matr);

    void transpose() const;

    void add_row_to_row(size_t src_num, size_t dst_num) const;

    void sub_row_to_row(size_t src_num, size_t dst_num) const;

    void sub_row_to_row(size_t src_num, size_t dst_num, T koef) const;

    void swap_rows(size_t row1_num, size_t row2_num) const;

    void mult_row_on_number(size_t row_num, T number) const;

    T operator()(size_t row_i, size_t column_i) const override;

    size_t size() const { return size_; }

    bool operator==(const Abstract_Matrix<T>& inp_rhs) const override;

    Square_Matrix& operator+=(const Square_Matrix& rhs);

    Square_Matrix& operator-=(const Square_Matrix<T>& rhs);

    Square_Matrix& operator*=(const Square_Matrix<T>& rhs);

    void print() const override;

    T determinant() const;
};


template <typename T>
class Matrix final : public Abstract_Matrix<T> {
private:
    size_t column_size_ = 0;
    size_t row_size_ = 0;
    T **data_ = nullptr;

    using Abstract_Matrix<T>::is_match;

    void new_data();

    static T** new_data_buf(size_t column_size, size_t row_size);

    void delete_data();
public:
    //matrix elems should be contained in inp_data row by row
    Matrix(const T *const *inp_data, size_t inp_column_size, size_t inp_row_size);

    Matrix(const std::vector<std::vector<T>>& input_rows);

    ~Matrix() override;

    Matrix(const Matrix& matr): Matrix(matr.data_, matr.column_size_, matr.row_size_) {}

    Matrix& operator= (const Matrix& matr)&;

    template<typename U> Matrix(const Matrix<U>& matr);

    Matrix(const Square_Matrix<T>& matr): Matrix(matr.data_, matr.size_, matr.size_) {}

    void transpose();

    void add_row_to_row(size_t src_num, size_t dst_num) const;

    void sub_row_to_row(size_t src_num, size_t dst_num) const;

    void sub_row_to_row(size_t src_num, size_t dst_num, T koef) const;

    void swap_rows(size_t row1_num, size_t row2_num) const;

    void mult_row_on_number(size_t row_num, T number) const;

    T operator()(size_t row_i, size_t column_i) const override;

    size_t row_size() const { return row_size_; }
    size_t column_size() const { return column_size_; }

    bool operator==(const Abstract_Matrix<T>& inp_rhs) const override;

    Matrix& operator+=(const Matrix& rhs);

    Matrix& operator-=(const Matrix<T>& rhs);

    Matrix& operator*=(const Matrix<T>& rhs);

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
void Square_Matrix<T>::new_data() {
    if(size_ > 0) {
        data_ = new T*[size_];
        for(size_t i = 0; i < size_; ++i) {
            data_[i] = new T[size_];
        }
    }
}

template<typename T>
T** Square_Matrix<T>::new_data_buf(size_t size) {
    if(size == 0) return nullptr;
    T** data = new T*[size];
    for(size_t i = 0; i < size; ++i) {
        data[i] = new T[size];
    }
    return data;
}

template<typename T>
void Square_Matrix<T>::delete_data() {
    if(size_ > 0) {
        for(size_t i = 0; i < size_; ++i) {
            delete [] data_[i];
        }
        delete [] data_;
    }
}



template<typename T>
Square_Matrix<T>::Square_Matrix(const T *const *inp_data, size_t inp_size):
    size_(inp_size)
{
    new_data();

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i][j] = inp_data[i][j];
        }
    }
}

template<typename T>
Square_Matrix<T>::Square_Matrix(const std::vector<std::vector<T>>& input_rows):
    size_(input_rows.size())
{
    for(const std::vector<T> &row : input_rows) {
        assert((row.size() == size_) && ("Invalid matrix size"));
    }

    new_data();

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i][j] = input_rows[i][j];
        }
    }
}

template<typename T>
Square_Matrix<T>::~Square_Matrix() {
    delete_data();
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator= (const Square_Matrix<T>& matr)& {
    if(this != &matr) {
        delete_data();

        size_ = matr.size();

        data_ = new_data();

        for(size_t i = 0; i < size_; ++i) {
            for(size_t j = 0; j < size_; ++j) {
                data_[i][j] = matr(i, j);
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

    new_data();

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i][j] = static_cast<T>(matr(i, j));
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
    assert((src_num < size_) && "invalid row number");
    assert((dst_num < size_) && "invalid row number");

    for(size_t i = 0; i < size_; ++i) {
        data_[dst_num][i] += data_[src_num][i];
    }
}

template<typename T>
void Square_Matrix<T>::sub_row_to_row(size_t src_num, size_t dst_num) const {
    assert((src_num < size_) && "invalid row number");
    assert((dst_num < size_) && "invalid row number");

    for(size_t i = 0; i < size_; ++i) {
        data_[dst_num][i] -= data_[src_num][i];
    }
}

template<typename T>
void Square_Matrix<T>::sub_row_to_row(size_t src_num, size_t dst_num, T koef) const {
    assert((src_num < size_) && "invalid row number");
    assert((dst_num < size_) && "invalid row number");

    for(size_t i = 0; i < size_; ++i) {
        data_[dst_num][i] -= (data_[src_num][i] * koef);
    }
}

template<typename T>
void Square_Matrix<T>::swap_rows(size_t row1_num, size_t row2_num) const {
    assert((row1_num < size_) && (row2_num < size_) && "invalid row_num");
    std::swap(data_[row1_num], data_[row2_num]);
}

template<typename T>
void Square_Matrix<T>::mult_row_on_number(size_t row_num, T number) const {
    assert((row_num < size_) && "invalid row_num");
    for(size_t i = 0; i < size_; ++i) {
        data_[row_num][i] *= number;
    }
}

template<typename T>
T Square_Matrix<T>::operator()(size_t row_i, size_t column_i) const {
    assert((row_i < size_) && "invalid column");
    assert((column_i < size_) && "invalid row");
    return data_[row_i][column_i];
}

template<typename T>
bool Square_Matrix<T>::operator==(const Abstract_Matrix<T>& inp_rhs) const {
    if(typeid(*this) != typeid(inp_rhs)) return false; //IDE?

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
Square_Matrix<T>& Square_Matrix<T>::operator+=(const Square_Matrix<T>& rhs) {
    assert((size_ == rhs.size()) && "different matrix sizes");
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i][j] += rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator-=(const Square_Matrix<T>& rhs) {
    assert((size_ == rhs.size()) && "different matrix sizes");
    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            data_[i][j] -= rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Square_Matrix<T>& Square_Matrix<T>::operator*=(const Square_Matrix<T>& rhs) { //naive algorithm (almost)
    assert((size_ == rhs.size()) && "different matrix sizes");
    if(size_ == 0) return *this;

    T** new_data = new_data_buf(size_);

    Square_Matrix new_rhs(rhs);
    new_rhs.transpose(); //for cache friendly

    for(size_t i = 0; i < size_; ++i) {
        for(size_t j = 0; j < size_; ++j) {
            new_data[i][j] = 0;
            for(size_t k = 0; k < size_; ++k) {
                new_data[i][j] += (data_[i][k] * new_rhs(j, k));
            }
        }
    }

    delete_data();
    data_ = new_data;

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



//----------------------------methods for Matrix----------------------------------------

template<typename T>
void Matrix<T>::new_data() {
    if((column_size_ > 0) && (row_size_ > 0)) {
        data_ = new T*[column_size_];
        for(size_t i = 0; i < column_size_; ++i) {
            data_[i] = new T[row_size_];
        }
    }
}

template<typename T>
T** Matrix<T>::new_data_buf(size_t column_size, size_t row_size) {
    if((column_size == 0) || (row_size == 0)) return nullptr;
    T** data = new T*[column_size];
    for(size_t i = 0; i < column_size; ++i) {
        data[i] = new T[row_size];
    }
    return data;
}

template<typename T>
void Matrix<T>::delete_data() {
    if((column_size_ > 0) && (row_size_ > 0)) {
        for(size_t i = 0; i < column_size_; ++i) {
            delete [] data_[i];
        }
        delete [] data_;
    }
}

template<typename T>
Matrix<T>::Matrix(const T *const *inp_data, size_t inp_column_size, size_t inp_row_size):
    column_size_(inp_column_size),
    row_size_(inp_row_size)
{
    new_data();

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            data_[i][j] = inp_data[i][j];
        }
    }
}

template<typename T>
Matrix<T>::Matrix(const std::vector<std::vector<T>>& input_rows):
    column_size_(input_rows.size()),
    row_size_((input_rows.size() > 0) ? input_rows[0].size() : 0)
{
    for(const std::vector<T> &row : input_rows) {
        assert((row.size() == row_size_) && ("Invalid matrix size"));
    }

    new_data();

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            data_[i][j] = input_rows[i][j];
        }
    }
}

template<typename T>
Matrix<T>::~Matrix() {
    delete_data();
}

template<typename T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& matr)& {
    if(this != &matr) {
        delete_data();

        column_size_ = matr.column_size();
        row_size_ = matr.row_size();

        new_data();

        for(size_t i = 0; i < column_size_; ++i) {
            for(size_t j = 0; j < row_size_; ++j) {
                data_[i][j] = matr(i, j);
            }
        }
    }

    return *this;
}

template<typename T>
template<typename U>
Matrix<T>::Matrix(const Matrix<U>& matr):
    column_size_(matr.column_size()),
    row_size_(matr.row_size())
{
    new_data();
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            data_[i][j] = static_cast<T>(matr(i, j));
        }
    }
}

template<typename T>
void Matrix<T>::transpose() {
    T** new_data = new_data_buf(row_size_, column_size_);

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            new_data[j][i] = data_[i][j];
        }
    }

    delete_data();

    std::swap(row_size_, column_size_);
    data_ = new_data;
}

template<typename T>
void Matrix<T>::add_row_to_row(size_t src_num, size_t dst_num) const {
    assert((src_num < column_size_) && "invalid row number");
    assert((dst_num < column_size_) && "invalid row number");

    for(size_t i = 0; i < row_size_; ++i) {
        data_[dst_num][i] += data_[src_num][i];
    }
}

template<typename T>
void Matrix<T>::sub_row_to_row(size_t src_num, size_t dst_num) const {
    assert((src_num < column_size_) && "invalid row number");
    assert((dst_num < column_size_) && "invalid row number");

    for(size_t i = 0; i < row_size_; ++i) {
        data_[dst_num][i] -= data_[src_num][i];
    }
}

template<typename T>
void Matrix<T>::sub_row_to_row(size_t src_num, size_t dst_num, T koef) const {
    assert((src_num < column_size_) && "invalid row number");
    assert((dst_num < column_size_) && "invalid row number");

    for(size_t i = 0; i < row_size_; ++i) {
        data_[dst_num][i] -= (data_[src_num][i] * koef);
    }
}

template<typename T>
void Matrix<T>::swap_rows(size_t row1_num, size_t row2_num) const {
    assert((row1_num < column_size_) && (row2_num < column_size_) && "invalid row_num");
    std::swap(data_[row1_num], data_[row2_num]);
}

template<typename T>
void Matrix<T>::mult_row_on_number(size_t row_num, T number) const {
    assert((row_num < column_size_) && "invalid row_num");
    for(size_t i = 0; i < row_size_; ++i) {
        data_[row_num][i] *= number;
    }
}

template<typename T>
T Matrix<T>::operator()(size_t column_i, size_t row_i) const {
    assert((row_i < row_size_) && "invalid column");
    assert((column_i < column_size_) && "invalid row");
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
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs) {
    assert((row_size_ == rhs.row_size()) && "different matrix sizes");
    assert((column_size_ == rhs.column_size()) && "different matrix sizes");
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            data_[i][j] += rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs) {
    assert((row_size_ == rhs.row_size()) && "different matrix sizes");
    assert((column_size_ == rhs.column_size()) && "different matrix sizes");
    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < row_size_; ++j) {
            data_[i][j] -= rhs(i, j);
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& rhs) { //naive algorithm (almost)
    assert((row_size_ == rhs.column_size()) && "invalid matrix sizes");

    if((row_size_ == 0) || (column_size_ == 0)) return *this;

    if(rhs.row_size() == 0) {
        delete_data();
        data_ = nullptr;
        row_size_ = 0;
        return *this;
    }

    T** new_data = new_data_buf(column_size_, rhs.row_size());

    Matrix new_rhs(rhs);
    new_rhs.transpose(); //for cache friendly

    for(size_t i = 0; i < column_size_; ++i) {
        for(size_t j = 0; j < new_rhs.column_size(); ++j) {
            new_data[i][j] = 0;
            for(size_t k = 0; k < row_size_; ++k) {
                new_data[i][j] += (data_[i][k] * new_rhs(j, k));
            }
        }
    }

    delete_data();
    row_size_ = rhs.row_size();
    data_ = new_data;

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
