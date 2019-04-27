#include <fstream>
#include <tsp_instance.hpp>
#include <two_matching_cut_gen.hpp>

int main() {
	std::ifstream in("../instances/pcb442.tsp");
	TSPInstance inst(in);
	std::vector<double> sol(inst.getEdgeCount());
	std::ifstream data("broken2matching.txt");
	std::string line;
	while (std::getline(data, line) && !line.empty()) {
		std::stringstream ss(line);
		size_t index;
		double val;
		ss >> index;
		ss.ignore();
		ss >> val;
		sol[index] = val;
	}
	int status;
	CPXENVptr env = CPXopenCPLEX(&status);
	if (status!=0) {
		throw std::runtime_error("Failed to open CPLEX environment: "+std::to_string(status));
	}
	TspLpData lpData(inst, nullptr);
	LinearProgram lp(env, "test", LinearProgram::minimize);
	lpData.setupBasicLP(lp);
	TwoMatchingCutGen cg(inst, lpData, false);
	std::cout << cg.validate(lp, sol, CutGenerator::valid) << std::endl;
	CPXcloseCPLEX(&env);
}