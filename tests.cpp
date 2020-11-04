#include <iostream>

#include "matrix.h"

using namespace matrix;

void test0() {
    Square_Matrix<int> m({{1, 1, 1},
                          {1, 1, 1},
                          {1, 1, 1}});

    int det = determinant(m);

    assert(det == 0);
    std::cout << "Test0 complete.\n";
}

void test1() {
    Square_Matrix<int> m({{2, 0, 1},
                          {1, 0, 1},
                          {1, 1, 1}});

    int det = determinant(m);

    assert(det == -1);
    std::cout << "Test1 complete.\n";
}
