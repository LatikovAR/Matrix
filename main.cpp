//This programm counts determinant
//Now it can work with types: double, int, long long int
//TEST - check all tests
//WORK - normal work with cin/cout

#define GENERATOR
#include <iostream>
#include <cstdlib>
#include <vector>

#include "matrix.h"

#ifdef TEST
#include "tests.h"
using namespace tests;
#endif

#ifdef GENERATOR
#include "generator.h"
using namespace generator;
#endif

int main()
{
#ifdef GENERATOR
    generate_matrix(100, 42);
    std::cout << "Generation complete\n";
#endif

#ifdef TEST
    unit_test0();
    unit_test1();
    unit_test2();
    unit_test3();
    unit_test4();
    unit_test5();
    unit_test6();
    unit_test7();

    test0();
    test1();
    test2();
    long_test(10);
#endif
#ifdef WORK
    size_t matr_size = 0;
    std::cin >> matr_size;
    int *data = new int[matr_size * matr_size]; //change type if nessesary
    for(size_t i = 0; i < matr_size * matr_size; ++i) {
        std::cin >> data[i];
    }

    matrix::Square_Matrix<int> matr(data, matr_size); //change type if nessesary
    auto det = matr.determinant();
    std::cout << det << std::endl;

    delete [] data;
#endif
    return 0;
}
