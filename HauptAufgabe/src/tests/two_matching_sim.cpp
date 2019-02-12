#include <fstream>
#include <tsp_instance.hpp>
#include <two_matching_cut_gen.hpp>

int main() {
	std::ifstream in("../instances/gr431.tsp");
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
	LinearProgram lp("test", LinearProgram::minimize);
	inst.setupBasicLP(lp);
	TwoMatchingCutGen cg(inst, false);
	std::cout << cg.validate(lp, sol) << std::endl;
}