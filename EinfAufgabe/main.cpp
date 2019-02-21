#include "graph.hpp"
#include "tspsolvers.hpp"
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <stdexcept>

int main(int argc, char** argv) {
	if (argc!=2) {
		std::cout << "Requires exactly one argument (input file name)" << std::endl;
		return 1;
	}
	std::ifstream input(argv[1]);
	if (!input) {
		std::cout << "File does not exist" << std::endl;
		return 2;
	}
	try {
		Graph g(input);
		Tour ref = tspsolvers::getReferenceTour(g);
		Tour greedy = tspsolvers::greedyTSP(g);
		std::cout << "Cost of reference tour: " << ref.length << std::endl;
		std::cout << "Cost of greedy tour: " << greedy.length << std::endl;
		greedy.print(std::cout);
	} catch (std::invalid_argument &x) {
		std::cout << "ERROR: " << x.what() << std::endl;
	}
}
