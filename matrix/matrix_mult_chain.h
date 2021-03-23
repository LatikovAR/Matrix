#pragma once

#include <iostream>
#include <vector>
#include <stdexcept>
#include <map>
#include <set>
#include <string>

#include "matrix.h"

namespace matrix {

template<typename T>
class Matrix_Chain final {
private:

    //Trace is a way of multiplication execution
    //contains numbers of "*" in some execution order
    class Trace final {
        std::vector<size_t> op_order_; //operations order
        size_t weight_; //number of trivial "*" operations

    public:
        Trace(): weight_(0) {}

        std::vector<size_t> op_order() const { return op_order_; }
        void write_op_order(std::vector<size_t> op_order) {
            op_order_ = std::move(op_order);
        }
        void push_op_num_back(size_t op_num) { op_order_.push_back(op_num); }

        size_t& weight() { return weight_; }
    };

    //inserting in vector isn't so fast as in list
    //but vector is better for elems access
    std::vector<Matrix<T>> chain_;

    //every bit in elem num performs used "*"
    //00000000 = 0 - nothing used
    //00000010 = 2 - used only second "*"
    //00000110 = 6 - used third and second "*"
    //11111111 - used all (answer here)
    std::vector<Trace> traces_;

    //num of using bits in traces numbers
    size_t bit_length_;

    //contains column sizes for all matrix in right order
    //and last matrix row size at the end;
    std::vector<size_t> matr_sizes;

    void make_first_trace_level();
    void make_trace_level(size_t i);
public:
    Matrix_Chain() = default;
    ~Matrix_Chain() = default;

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
    make_first_trace_level(); //only resize traces vector

    for(size_t i = 0; i < traces_.size() - 1; ++i) {
        make_trace_level(i); //build all possible ways from this condition
    }

    auto last_elem = traces_.rbegin();
    return last_elem->op_order();
}

template<typename T>
void Matrix_Chain<T>::make_first_trace_level() {
    traces_.clear();

    size_t c_size = chain_.size() - 1;
    bit_length_ = c_size;

    //calculating 2^bit_length
    size_t t_num = 1;
    size_t multiplicator = 2;
    while(c_size > 0) {
        if(c_size % 2 == 0) {
            multiplicator *= multiplicator;
            c_size /= 2;
        }
        else {
            t_num *= multiplicator;
            c_size--;
        }
    }

    traces_.resize(t_num);

    matr_sizes.resize(chain_.size() + 1);
    for(size_t i = 0; i < chain_.size(); ++i) {
        matr_sizes[i] = chain_[i].column_size();
    }
    matr_sizes[chain_.size()] = chain_[chain_.size() - 1].row_size();
}

template<typename T>
void Matrix_Chain<T>::make_trace_level(size_t cur_pos) {
    size_t bot_size, mid_size, top_size;
    size_t mid_num;
    bool mid_found = false;

    bot_size = matr_sizes[0];
    for(size_t mult_num = 0; mult_num < bit_length_; ++mult_num) {
        if((cur_pos >> mult_num) % 2 == 0) {

            if(mid_found == false) {
                mid_found = true;
                mid_num = mult_num;
                mid_size = matr_sizes[mult_num + 1];
            }

            else {
                top_size = matr_sizes[mult_num + 1];

                size_t next_pos = cur_pos | (1 << mid_num);
                size_t weight_inc = bot_size * mid_size * top_size;

                if((traces_[next_pos].weight() == 0) ||
                   (traces_[next_pos].weight() > traces_[cur_pos].weight() + weight_inc))
                {
                    traces_[next_pos].weight() = traces_[cur_pos].weight() + weight_inc;
                    traces_[next_pos].write_op_order(std::move(traces_[cur_pos].op_order()));
                    traces_[next_pos].push_op_num_back(mid_num);
                }

                bot_size = mid_size;
                mid_size = top_size;
                mid_num = mult_num;
            }
        }
    }

    if(mid_found == true) {
        top_size = *matr_sizes.rbegin();

        size_t next_pos = cur_pos | (1 << mid_num);
        size_t weight_inc = bot_size * mid_size * top_size;

        if((traces_[next_pos].weight() == 0) ||
           (traces_[next_pos].weight() > traces_[cur_pos].weight() + weight_inc))
        {
            traces_[next_pos].weight() = traces_[cur_pos].weight() + weight_inc;
            traces_[next_pos].write_op_order(std::move(traces_[cur_pos].op_order()));
            traces_[next_pos].push_op_num_back(mid_num);
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
