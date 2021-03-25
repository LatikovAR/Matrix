cmake_minimum_required(VERSION 3.5)

project(matrix_mult LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(${CMAKE_SOURCE_DIR}/matrix)

include_directories(${CMAKE_SOURCE_DIR}/mem_storage ${CMAKE_SOURCE_DIR}/dynamic_array)

add_executable(matrix_mult ${CMAKE_SOURCE_DIR}/main.cpp)

target_link_libraries(matrix_mult ${CMAKE_SOURCE_DIR}/matrix)
