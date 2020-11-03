#include <cassert>
#include <vector>
#include <iostream>

#include "tests.h"
#include "matrix.h"

using namespace matrix;

void unit_test0() {
    Matrix<int> m({{2, 0},
                   {1, 5},
                   {-1, -2}});

    Matrix<int> m_new = m.transpose();
    Matrix<int> m_new_sample({{2, 1, -1},
                              {0, 5, -2}});

    assert(m_new == m_new_sample);
    std::cout << "Unit_test0 complete.\n";
}

void unit_test1() {
    Matrix<double> m({{2.1, 0.5},
                   {1.8, 5.4},
                   {-1.9, -2.6}});

    Matrix<double> m_new = m.transpose();
    Matrix<double> m_new_sample({{2.1, 1.8, -1.9},
                              {0.5, 5.4, -2.6}});

    assert(m_new == m_new_sample);
    std::cout << "Unit_test1 complete.\n";
}
