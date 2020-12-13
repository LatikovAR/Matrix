project(circuit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCE_EXE main.cpp circuit.cpp input_output/input_output.cpp)

set(SOURCE_LIB matrix/matrix.cpp)

add_library(matrix STATIC ${SOURCE_LIB})

add_executable(circuit_solver ${SOURCE_EXE})

target_link_libraries(circuit_solver matrix)	
