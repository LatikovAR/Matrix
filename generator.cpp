#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>

#include "matrix.h"

using namespace matrix;

namespace {

//it generates L and U matrices of LU-decomposition, multiply and shuffle
Square_Matrix<long long> gen_matrix(size_t size, long long det) {
    long long *data = new long long[size * size];
    srand(static_cast<unsigned>(clock()));

    //gen L
    for(size_t i = 0; i < size; ++i) {
        for(size_t j = 0; j < size; ++j) {
            long long buf = 0;
            if(i > j) {
                buf = rand() % 2;
            }
            if(i == j) {
                buf = 1;
            }
            if(i < j) {
                buf = 0;
            }
            data[i * size + j] = buf;
        }
    }
    Square_Matrix<long long> L(data, size);

    //gen U
    for(size_t i = 0; i < size; ++i) {
        for(size_t j = 0; j < size; ++j) {
            long long buf = 0;
            if(i < j) {
                buf = rand() % 2;
            }
            if(i == j) {
                if(i == 0) {
                    buf = det; //not so interesting, but easy :)
                }
                else {
                    buf = 1;
                }
            }
            if(i > j) {
                buf = 0;
            }
            data[i * size + j] = buf;
        }
    }
    Square_Matrix<long long> U(data, size);

    //multiply
    Square_Matrix<long long> res = L * U;

    //shuffle (make matrix more "interesting"
    if(size > 1) {
        for(size_t i = 0; i < (size / 10 + 1); ++i) {
            size_t src = static_cast<size_t>(abs(rand())) % size;
            size_t dst;
            do {
                dst = static_cast<size_t>(abs(rand())) % size;
            } while(src == dst);
            res.sub_row_to_row(src, dst);
        }
    }

    return res;
}

//better do it on templates but this way easier
void print_matrix_in_file(char* filename, Square_Matrix<long long> m) {
    FILE *f = fopen(filename, "w");
    fprintf(f, "%llu ", static_cast<unsigned long long>(m.size()));
    for(size_t i = 0; i < m.size(); ++i) {
        for(size_t j = 0; j < m.size(); ++j) {
            fprintf(f, "%lld ", m(i, j));
        }
    }
    fclose(f);
}

}

namespace generator {
void generate_matrix(size_t size, long long det) {
    Square_Matrix<long long> m = gen_matrix(size, det);

    std::cout << m.determinant() << std::endl;
    char filename[100] = "001.dat";
    print_matrix_in_file(filename, m);

    FILE *f = fopen("001.ans", "w");
    fprintf(f, "%lld", det);
    fclose(f);
}
}
