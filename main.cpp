#include <iostream>

#include "circuit.h"

int main() {
    std::vector<circuit::Edge_Info> edges_info = circuit::Edge_Info::input_edges_info();
    circuit::Circuit circuit(edges_info);

    if(circuit.validity() == true) {
        circuit.print_circuit();
    }
    else {
        std::cout << "Input circuit can't be solved.\n";
        std::cout << "Maybe you set some R <= 0.0.\n";
    }
    return 0;
}
