//Matrix_Chain implements searching the right order of matrix multiplication.
//Imagine that we have matrices multiplication A * B
//A sizes is (n x m), B sizes is (m x l)
//so we have n * m * l trivial "*" operations
//A * B * C = (A * B) * C = A * (B * C)
//so it's good to choose variant with minimal trivial "*" number.
//
//Introduce notations:
//weight - number of trivial "*"
//A(i, j) - result of multiplication all matrices from i to j number
//
//Obviously that A(i, i) - it is A with number i without any multiplication
//A(i, j) = A(i, k) * A(k + 1, j) where i <= k < j
//so used method of dynamic programming for all A(i, j), i <= j

#pragma once

#include <iostream>
#include <vector>
#include <stdexcept>
#include <map>
#include <set>
#include <string>

#include "matrix.h"

namespace matrix {

namespace {

//Trace is a way of multiplication execution
//contains numbers of "*" in some execution order
class Trace final {
    std::vector<size_t> op_order_; //operations order

public:
    size_t weight = 0; //number of trivial "*" operations

    std::vector<size_t> op_order() const { return op_order_; }

    //places first_part before second_part and new_op in the end
    void build_op_order(std::vector<size_t> first_part,
                        std::vector<size_t> second_part,
                        size_t new_op)
    {
        op_order_.clear();
        op_order_.reserve(first_part.size() + second_part.size() + 1);

        op_order_.insert(op_order_.end(), first_part.begin(), first_part.end());
        op_order_.insert(op_order_.end(), second_part.begin(), second_part.end());
        op_order_.insert(op_order_.end(), new_op);
    }
};

} //namespace

template<typename T>
class Matrix_Chain final {
private:
    //inserting in vector isn't so fast as in list
    //but vector is better for elems access
    std::vector<Matrix<T>> chain_;

    //traces is conveniently to store in matrix
    //A(i, j) => i and j - indices in matrix
    //Symmetric_Matrix because i <= j
    Symmetric_Matrix<Trace> traces_;

    //contains column sizes for all matrix in right order
    //and last matrix row size at the end;
    std::vector<size_t> matr_sizes;

    void make_start_trace_conditions();
    void make_trace_condition(size_t i, size_t j);
public:

    size_t matrix_num() const { return chain_.size(); }
    const Matrix<T>& ret_matrix(size_t number) const;

    //insert matrix in chain before matrix with number "pos"
    //so added matrix have number "pos" after inserting
    //allowed pos == chain_.size() for inserting at the end
    void insert_matrix(Matrix<T> matr, size_t pos);

    std::vector<size_t> compute_optimal_trace();

    static void print_optimal_trace_with_numbers(const std::vector<size_t>&);
    static void print_optimal_trace_with_brackets(const std::vector<size_t>&);

    Matrix<T> mult_chain_optimal(const std::vector<size_t>& order) const;
    Matrix<T> mult_chain_naive() const;
};

//for using matrix
template<>
bool Symmetric_Matrix<Trace>::operator==(const Abstract_Matrix<Trace>& inp_rhs) const {
    return false;
}

template<>
void Symmetric_Matrix<Trace>::print() const {}

template<typename T>
const Matrix<T>& Matrix_Chain<T>::ret_matrix(size_t number) const {
    if(number >= chain_.size())
        throw std::invalid_argument("Matrix_chain: no matrix with this number");

    return chain_[number];
}

template<typename T>
void Matrix_Chain<T>::insert_matrix(Matrix<T> matr, size_t pos) {
    if(pos > chain_.size())
        throw std::invalid_argument("Matrix_chain: can't be added matrix at this pos");

    if(pos != 0) {
        if(chain_[pos - 1].row_size() != matr.column_size()) {
            throw std::invalid_argument("Matrix_chain: invalid inserting matrix size");
        }
    }

    if(pos != chain_.size()) {
        if(chain_[pos].column_size() != matr.row_size()) {
            throw std::invalid_argument("Matrix_chain: invalid inserting matrix size");
        }
    }

    chain_.insert(chain_.begin() + pos, std::move(matr));
}

template<typename T>
std::vector<size_t> Matrix_Chain<T>::compute_optimal_trace() {
    make_start_trace_conditions(); //build A(i, i) conditions - without "*"

    for(size_t i = 1; i < traces_.size(); ++i) {
        for(size_t j = 0; j + i < traces_.size(); ++j) {
            make_trace_condition(j, j + i); //check all possible ways to this condition and choose the best
        }
    }

    return traces_(0, traces_.size() - 1).op_order();
}

template<typename T>
void Matrix_Chain<T>::make_start_trace_conditions() {
    traces_ = Symmetric_Matrix<Trace>(chain_.size());

    matr_sizes.resize(chain_.size() + 1);
    for(size_t i = 0; i < chain_.size(); ++i) {
        matr_sizes[i] = chain_[i].column_size();
    }
    matr_sizes[chain_.size()] = chain_[chain_.size() - 1].row_size();
}

template<typename T>
void Matrix_Chain<T>::make_trace_condition(size_t i, size_t j) {
    for(size_t k = i; k < j; ++k) {
        size_t new_weight = traces_(i, k).weight +
                            matr_sizes[i] * matr_sizes[k + 1] * matr_sizes[j + 1] +
                            traces_(k + 1, j).weight;

        if((traces_(i, j).weight > new_weight) || (traces_(i, j).weight == 0)) {
            traces_(i, j).weight = new_weight;
            traces_(i, j).build_op_order(std::move(traces_(i, k).op_order()),
                                         std::move(traces_(k + 1, j).op_order()),
                                         k);
        }
    }
}

template<typename T>
void Matrix_Chain<T>::print_optimal_trace_with_numbers(const std::vector<size_t>& trace) {
    for(const auto& elem : trace) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

template<typename T>
void Matrix_Chain<T>::print_optimal_trace_with_brackets(const std::vector<size_t>& trace) {
    std::map<size_t, std::string> output;

    for(size_t i = 0; i < trace.size() + 1; ++i) {
        output.emplace(i, "M");
    }

    for(const auto& op_num : trace) {
        auto it1 = output.find(op_num);
        assert(it1 != output.end());
        auto it2 = it1;
        it2++;
        assert(it2 != output.end());
        std::string& str1 = it1->second;
        std::string& str2 = it2->second;

        bool lhs_brace = (str1 != "M");
        bool rhs_brace = (str2 != "M");

        if(lhs_brace) str1 = "(" + str1 + ")";
        if(rhs_brace) str2 = "(" + str2 + ")";

        str2 = str1 + " * " + str2;

        output.erase(it1);
    }

    assert(output.size() == 1);
    auto it = output.find(trace.size());

    std::cout << it->second << std::endl;
}

template<typename T>
Matrix<T> Matrix_Chain<T>::mult_chain_optimal(const std::vector<size_t>& order) const {
    if(chain_.size() == 0)
        throw std::runtime_error("Mult chain: no matrices");

    if(order.size() + 1 != chain_.size())
        throw std::runtime_error("Mult chain: not suitable order param");

    std::map<size_t, Matrix<T>> expr;
    for(size_t i = 0; i < chain_.size(); ++i) {
        expr.emplace(i, chain_[i]);
    }

    for(const auto& mult_num : order) {

        auto it1 = expr.find(mult_num);
        if(it1 == expr.end())
            throw std::runtime_error("Mult chain: not suitable order param");
        auto it2 = it1;
        it2++;
        if(it1 == expr.end())
            throw std::runtime_error("Mult chain: not suitable order param");

        it2->second = it1->second * it2->second;
        expr.erase(it1);
    }

    assert(expr.size() == 1);
    auto it = expr.find(chain_.size() - 1);
    assert(expr.end() != it);
    return it->second;
}

template<typename T>
Matrix<T> Matrix_Chain<T>::mult_chain_naive() const {
    if(chain_.size() == 0)
        throw std::runtime_error("Mult chain: no matrices");

    Matrix<T> buf{chain_[0]};
    for(size_t i = 1; i < chain_.size(); ++i) {
        buf = buf * chain_[i];
    }

    return buf;
}

} //namespace matrix
