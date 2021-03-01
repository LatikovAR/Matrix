cmake_minimum_required(VERSION 3.5)

project(matrix_mult LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(Matrix)

include_directories(mem_storage dynamic_array)

add_executable(matrix_mult main.cpp matrix)
