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

    m.transpose();
    Matrix<int> m_new({{2, 1, -1},
                       {0, 5, -2}});

    assert(m == m_new);
    std::cout << "Unit_test0 complete.\n";
}

void unit_test1() {
    Matrix<double> m({{2.1, 0.5},
                   {1.8, 5.4},
                   {-1.9, -2.6}});

    m.transpose();
    Matrix<double> m_new({{2.1, 1.8, -1.9},
                          {0.5, 5.4, -2.6}});

    assert(m_new == m);
    std::cout << "Unit_test1 complete.\n";
}

void unit_test2() {
    Square_Matrix<double> m({{2.1, 0.5},
                             {1.8, 5.4}});

    m.transpose();
    Square_Matrix<double> m_new({{2.1, 1.8},
                                 {0.5, 5.4}});

    assert(m_new == m);
    std::cout << "Unit_test2 complete.\n";
}

void unit_test3() {
    Square_Matrix<double> m({});

    m.transpose();
    Square_Matrix<double> m_new({});

    assert(m_new == m);
    std::cout << "Unit_test3 complete.\n";
}

void unit_test4() {
    Square_Matrix<int> m({{2, 0},
                          {1, 5}});

    assert(m(0, 0) == 2);
    std::cout << "Unit_test4 complete.\n";
}

void unit_test5() {
    Matrix<int> m({{2, 0},
                   {1, 5},
                   {-1, -2}});

    m.add_row_to_row(1, 0);
    m.add_row_to_row(2, 1);
    Matrix<int> m_new({{3, 5},
                       {0, 3},
                       {-1, -2}});

    assert(m == m_new);
    std::cout << "Unit_test5 complete.\n";
}

void unit_test6() {
    Square_Matrix<double> m({{1.0, 0.0},
                             {0.0, 1.0}});

    double det = determinant(m);

    assert(det == 1.0);
    std::cout << "Unit_test6 complete.\n";
}

void unit_test7() {
    Square_Matrix<int> m({{1, 1},
                          {1, 1}});

    int det = determinant(m);

    assert(det == 0);
    std::cout << "Unit_test7 complete.\n";
}
