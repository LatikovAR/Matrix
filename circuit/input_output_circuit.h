#pragma once

#include "../matrix/matrix.h"

using namespace matrix;

namespace input_output {

class Input_Data final {
private:
    Symmetric_Matrix<bool> topology_; //true if edge exists
    Symmetric_Matrix<double> R_;
    Symmetric_Matrix<double> U_;
public:
    Input_Data();
    const Symmetric_Matrix<bool>& topology()& { return topology_; }
    const Symmetric_Matrix<double>& R()& { return R_; }
    const Symmetric_Matrix<double>& U()& { return U_; }
    void print_solution(const Symmetric_Matrix<double>& I);
};
}
