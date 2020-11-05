//This programm counts determinant
//Now it can work with types: double, int, long long int
//TEST - check all tests
//WORK - normal work with cin/cout

#define TEST
#include <iostream>
#include <cstdlib>
#include <vector>

#include "matrix.h"
#include "tests.h"

int main()
{
#ifdef TEST
    unit_test0();
    unit_test1();
    unit_test2();
    unit_test3();
    unit_test4();

    test0();
    test1();
    test2();
    long_test(10);
#endif
#ifdef WORK
    size_t matr_size = 0;
    std::cin >> matr_size;
    std::vector<std::vector<int>> rows; //change type if nessesary
    rows.resize(matr_size);
    for(size_t i = 0; i < matr_size; ++i) {
        rows[i].resize(matr_size);
    }
    for(size_t i = 0; i < matr_size; ++i) {
        for(size_t j = 0; j < matr_size; ++j) {
            std::cin >> rows[i][j];
        }
    }

    matrix::Square_Matrix<int> matr(rows); //change type if nessesary
    int det = matrix::determinant(matr);
    std::cout << det << std::endl;
#endif
    return 0;
}
