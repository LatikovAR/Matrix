#include <cassert>
#include <vector>
#include <iostream>

#include "tests.h"
#include "matrix.h"

using namespace matrix;

namespace tests {

void unit_test0() {
    Square_Matrix<double> m({{2.1, 0.5},
                             {1.8, 5.4}});

    m.transpose();
    Square_Matrix<double> m_new({{2.1, 1.8},
                                 {0.5, 5.4}});

    assert(m_new == m);
    std::cout << "Unit_test0 complete.\n";
}

void unit_test1() {
    Square_Matrix<double> m({});

    m.transpose();
    Square_Matrix<double> m_new({});

    assert(m_new == m);
    std::cout << "Unit_test1 complete.\n";
}

void unit_test2() {
    Square_Matrix<int> m({{2, 0},
                          {1, 5}});

    assert(m(0, 0) == 2);
    std::cout << "Unit_test2 complete.\n";
}

void unit_test3() {
    Square_Matrix<double> m({{1.0, 0.0},
                             {0.0, 1.0}});

    double det = m.determinant();

    assert(det == 1.0);
    std::cout << "Unit_test3 complete.\n";
}

void unit_test4() {
    Square_Matrix<int> m({{1, 1},
                          {1, 1}});

    int det = m.determinant();

    assert(det == 0);
    std::cout << "Unit_test4 complete.\n";
}

void unit_test5() {
    Square_Matrix<int> m1({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m2({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m3({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m1_res({{3, 3},
                               {3, 3}});

    Square_Matrix<int> m2_res({{2, 2},
                               {2, 2}});

    Square_Matrix<int> m3_res({{1, 1},
                               {1, 1}});

    m1 += m2 += m3;

    assert(m3 == m3_res);
    assert(m2 == m2_res);
    assert(m1 == m1_res);
    std::cout << "Unit_test5 complete.\n";
}

void unit_test6() {
    Square_Matrix<int> m1({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m2({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m3({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m1_res({{-1, -1},
                               {-1, -1}});

    Square_Matrix<int> m2_res({{2, 2},
                               {2, 2}});

    Square_Matrix<int> m3_res({{1, 1},
                               {1, 1}});

    m1 -= m2 *= m3;

    assert(m3 == m3_res);
    assert(m2 == m2_res);
    assert(m1 == m1_res);
    std::cout << "Unit_test6 complete.\n";
}

void unit_test7() {
    Square_Matrix<int> m1({{1, 0},
                           {0, 1}});

    Square_Matrix<int> m2({{1, 1},
                           {1, 1}});

    Square_Matrix<int> m3 = m1 * m2;

    Square_Matrix<int> m1_res({{1, 0},
                               {0, 1}});

    Square_Matrix<int> m2_res({{1, 1},
                               {1, 1}});

    Square_Matrix<int> m3_res({{1, 1},
                               {1, 1}});

    assert(m3 == m3_res);
    assert(m2 == m2_res);
    assert(m1 == m1_res);
    std::cout << "Unit_test7 complete.\n";
}
}
