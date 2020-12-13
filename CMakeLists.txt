project(circuit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(input_output matrix mem_storage dynamic_array)

set(SOURCE_EXE main.cpp circuit.cpp input_output/input_output.cpp matrix/matrix.cpp)

add_executable(circuit_solver ${SOURCE_EXE})
