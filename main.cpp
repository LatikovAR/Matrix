#include <iostream>

#include "circuit.h"
#include "input_output/input_output.h"
//#include "matrix/matrix_test.h"

int main() {
    //matrix_tests::matrix_test();
    input_output::Input_Data input_data;
    try {
        circuit::Circuit circuit(input_data.topology(),
                                 input_data.R(), input_data.U());
        input_data.print_solution(circuit.solve());
    }
    catch (circuit::Circuit::invalid_circuit& exc) {
        std::cout << exc.what() << std::endl;
    }

    return 0;
}
