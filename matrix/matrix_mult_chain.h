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

        //-------------------TRACE_MEMORY----------------
        std::vector<size_t> op_order_; //operations order
        size_t weight_; //number of trivial "*" operations

        //----------------FOR_CONTINUE_TRACE-------------
        //not used operation nums in right order like keys
        //matrix sizes around key operation like second param
        //-1 key - first matrix column size
        //last key - last matrix row size
        std::map<int, size_t> unused_ops;

    public:
        Trace(const std::vector<Matrix<T>>& chain): weight_(0) {
            size_t op_num = chain.size() - 1;

            op_order_.reserve(op_num);

            unused_ops.emplace(-1, chain[0].column_size());
            for(size_t i = 0; i < op_num; ++i) {
                unused_ops.emplace(i, chain[i].row_size());
            }
            unused_ops.emplace(op_num, chain[op_num].row_size());
        }

        std::vector<size_t> op_order() const { return op_order_; }

        size_t weight() const { return weight_; }

        void use_op(size_t op_num) {
            const auto& op = unused_ops.find(op_num);

            //no this op
            if(op == unused_ops.end())
                throw std::runtime_error("Trace: can't use operation with this num");
            auto prev_op = op;
            prev_op--;
            auto next_op = op;
            next_op++;

            op_order_.push_back(op_num);

            weight_ += (prev_op->second * op->second * next_op->second);

            unused_ops.erase(op);
        }

        size_t predict_weight(size_t op_num) const {
            const auto& op = unused_ops.find(op_num);

            //no this op
            if(op == unused_ops.end())
                throw std::runtime_error("Trace: can't use operation with this num");
            auto prev_op = op;
            prev_op--;
            auto next_op = op;
            next_op++;

            return weight_ + (prev_op->second * op->second * next_op->second);
        }
    };

    //inserting in vector isn't so fast as in list
    //but vector is better for elems access
    std::vector<Matrix<T>> chain_;

    //trace_levels_ - vector with calculation traces of multiplication
    //every level contains traces with the same length
    //first level <=> length = 0
    //std::map<std::set<size_t>, Trace> - map from one level traces
    //Key = std::set<size_t> - consists of calculated "*" numbers
    //for every Key stores the shortest trace
    std::vector<std::map<std::set<size_t>, Trace>> traces_levels_;

    std::map<std::set<size_t>, Trace> optimal_trace_;

    void make_first_trace_level();
    void make_trace_level(std::map<std::set<size_t>, Trace>& new_level,
                          const std::map<std::set<size_t>, Trace>& prev_level);
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
    make_first_trace_level();

    for(size_t i = 1; i < chain_.size(); ++i) {
        make_trace_level(traces_levels_[i], traces_levels_[i - 1]);
    }

    const auto& last_level = *traces_levels_.rbegin();
    assert(last_level.size() == 1);
    const Trace& final_trace = (*(last_level.begin())).second;
    return final_trace.op_order();
}

template<typename T>
void Matrix_Chain<T>::make_first_trace_level() {
    traces_levels_.clear();
    traces_levels_.resize(chain_.size());

    auto& first_level = traces_levels_[0];
    std::set<size_t> s;
    Trace tr(chain_); //other traces will be built over this
    first_level.emplace(std::move(s), std::move(tr));
}

template<typename T>
void Matrix_Chain<T>::make_trace_level(std::map<std::set<size_t>, Trace>& new_level,
                                       const std::map<std::set<size_t>, Trace>& prev_level)
{
    new_level.clear();

    for(const auto& prev_elem : prev_level) {
        const std::set<size_t>& prev_cond = prev_elem.first;
        const Trace& prev_trace = prev_elem.second;

        for(size_t mult_num = 0; mult_num < chain_.size() - 1; ++mult_num) {
            if(prev_cond.count(mult_num) == 0) { //if this "*" wasn't used before
                std::set<size_t> new_cond{prev_cond};
                new_cond.insert(mult_num);

                if(new_level.count(new_cond) == 0) {
                    Trace new_trace{prev_trace};
                    new_trace.use_op(mult_num);
                    new_level.emplace(std::move(new_cond), std::move(new_trace));
                }

                else {
                    assert(new_level.count(new_cond) == 1);

                    Trace& op_trace = new_level.find(new_cond)->second;
                    size_t new_weight = prev_trace.predict_weight(mult_num);

                    if(op_trace.weight() > new_weight) {
                        Trace new_trace{prev_trace};
                        new_trace.use_op(mult_num);
                        op_trace = std::move(new_trace);
                    }
                }
            }
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
    //not best algorithm but no matter here
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

} //namespace matrix
