cmake_minimum_required(VERSION 3.5)

project(libmatrix LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(matrix STATIC
    matrix.cpp
    matrix_mult_chain.cpp
    matrix_chain_test.cpp)
