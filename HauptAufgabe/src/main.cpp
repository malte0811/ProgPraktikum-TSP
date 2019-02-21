#include <tsp_solvers.hpp>
#include <tsp_instance.hpp>
#include <fstream>

int main(int argc, char** argv) {
	if (argc!=2 && argc!=3) {
		std::cout << "Arguments: <input file name> [<output file name>]" << std::endl;
		return 1;
	}
	std::ifstream in(argv[1]);
	if (!in) {
		std::cout << "File does not exist: " << argv[1] << std::endl;
		return 1;
	}
	try {
		int status;
		CPXENVptr env = CPXopenCPLEX(&status);
		if (status!=0) {
			throw std::runtime_error("Failed to open CPLEX environment: "+std::to_string(status));
		}
		TSPInstance inst(in);
		in.close();
		TSPSolution initial = tspsolvers::solveGreedy(inst);
		TSPSolution optimal = tspsolvers::solveLP(inst, &initial, env);
		if (argc==2) {
			optimal.write(std::cout);
		} else {
			std::ofstream out(argv[2]);
			if (!out) {
				std::cout << "Could not create/write to output file: " << argv[2] << std::endl;
				return 1;
			}
			optimal.write(out);
			out.close();
		}
		CPXcloseCPLEX(&env);
	} catch (std::runtime_error& err) {
		std::cout << "Error: " << err.what() << std::endl;
	}
}