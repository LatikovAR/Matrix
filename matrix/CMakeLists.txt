cmake_minimum_required(VERSION 3.5)

project(Lib_Matrix LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(matrix STATIC matrix.cpp matrix_mult_chain.cpp)
