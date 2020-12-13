#include <list>
#include <sstream>
#include <string>
#include <cstdlib>

#include "input_output.h"
#include "matrix/matrix.h"

namespace  {
struct Edge_Info {
    size_t begin_;
    size_t end_;
    double R_;
    double U_;

    Edge_Info(size_t begin, size_t end, double R, double U):
        begin_(begin), end_(end), R_(R), U_(U) {}

    //for correct input of the edges
    static bool is_prev_input_ok;
    static bool is_eof;

    static Edge_Info input_edge_info();
    static std::list<Edge_Info> input_edges_info();
    static size_t max_vert_num(const std::list<Edge_Info>& edges);
};

bool Edge_Info::is_prev_input_ok = true;
bool Edge_Info::is_eof = false;

//miss num symbols in the stream
//ignores space symbols
//0 - ok, -1 - wrong symbol or other problem
int miss_symbols(std::stringstream& stream, char symbol, int num = 1) {
    while(num > 0) {
        char rec_sym;
        stream >> rec_sym;
        if(isspace(rec_sym)) continue;

        if(!(stream.good()) || (symbol != rec_sym)) {
            return -1;
        }

        num--;
    }

    return 0;
}

Edge_Info Edge_Info::input_edge_info() {
    size_t begin = 0;
    size_t end = 0;
    double R = 0;
    double U = 0;

    std::string buf = {};

    is_prev_input_ok = true;
    is_eof = false;

    std::getline(std::cin, buf, '\n');
    if(std::cin.eof()) {
        is_eof = true;
        if(buf == "") { //sometimes it can be useful
            is_prev_input_ok = false;
            return Edge_Info(0, 0, 0.0, 0.0);
        }
    }
    if(std::cin.fail() && !(std::cin.eof())) {
        std::cout << "Warning: invalid str input format\n";
        is_prev_input_ok = false;
        return Edge_Info(0, 0, 0.0, 0.0);
    }

    std::stringstream s_buf;

    s_buf << buf;

    s_buf >> begin;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(0, 0, 0.0, 0.0);
    }

    if(miss_symbols(s_buf, '-', 2) < 0) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, 0, 0.0, 0.0);
    }

    s_buf >> end;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, 0, 0.0, 0.0);
    }

    if(miss_symbols(s_buf, ',') < 0) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, end, 0.0, 0.0);
    }

    s_buf >> R;
    if(!(s_buf.good())) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, end, 0.0, 0.0);
    }

    if(miss_symbols(s_buf, ';') < 0) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, end, 0.0, 0.0);
    }


    s_buf >> U;
    if(!(s_buf.good())) {
        if(s_buf.eof()) return Edge_Info(begin, end, R, 0.0); //correct input
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, end, R, 0.0);
    }

    if(miss_symbols(s_buf, 'V') < 0) {
        std::cout << "Warning: invalid input format\n";
        is_prev_input_ok = false;
        return Edge_Info(begin, end, R, 0.0);
    }

    char trash;
    s_buf >> trash;
    if(s_buf.eof()) return Edge_Info(begin, end, R, U); //correct input

    std::cout << "Warning: invalid input format\n";
    is_prev_input_ok = false;
    return Edge_Info(begin, end, R, U);
}

std::list<Edge_Info> Edge_Info::input_edges_info() {
    std::list<Edge_Info> edges;

    while(Edge_Info::is_eof == false) {
        Edge_Info edge_info = Edge_Info::input_edge_info();
        if(Edge_Info::is_prev_input_ok) {
            edges.push_back(edge_info);
        }
    }

    return edges;
}

size_t Edge_Info::max_vert_num(const std::list<Edge_Info>& edges) {
    size_t max_edge_num = 0;
    for(const Edge_Info edge : edges) {
        if(edge.begin_ > max_edge_num) max_edge_num = edge.begin_;
        if(edge.end_ > max_edge_num) max_edge_num = edge.end_;
    }
    return max_edge_num;
}

} //ending namespace

namespace input_output {

Input_Data::Input_Data() {
    std::list<Edge_Info> edges = Edge_Info::input_edges_info();

    size_t matr_size = Edge_Info::max_vert_num(edges) + 1;
    topology_ = Symmetric_Matrix<bool>(matr_size);
    R_ = Symmetric_Matrix<double>(matr_size);
    U_ = Symmetric_Matrix<double>(matr_size);

    for(size_t i = 0; i < matr_size; ++i) {
        for(size_t j = i; j < matr_size; ++j) {
            topology_(i, j) = false;
        }
    }

    for(const Edge_Info edge : edges) {
        topology_(edge.begin_, edge.end_) = true;
        R_(edge.begin_, edge.end_) = edge.R_;
        U_(edge.begin_, edge.end_) = edge.U_;
        if(edge.begin_ > edge.end_)
            U_(edge.begin_, edge.end_) *= -1.0;
    }
}

void Input_Data::print_solution(const Symmetric_Matrix<double>& I) {
    for(size_t i = 0; i < topology_.size(); ++i) {
        for(size_t j = i; j < topology_.size(); ++j) {
            if(topology_(i, j)) {
                std::cout << i << " -- " << j << " ";
                std::cout << I(i, j) << " A\n";
            }
        }
    }
}
} //ending namespace input_output
