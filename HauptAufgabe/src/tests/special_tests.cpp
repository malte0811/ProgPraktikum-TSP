#include <string>
#include <fstream>
#include <tsp_instance.hpp>
#include <subtour_cut_gen.hpp>
#include <two_matching_cut_gen.hpp>
#include <branch_and_cut.hpp>
#include <ctime>
#include <tsp_solution.hpp>
#include <tsp_solvers.hpp>

int main() {
	const std::vector<std::pair<std::string, cost_t>> instances{
			{"test_3",        3},
			{"test_4",        4},
			{"test_5",        5},
			{"large_weights", 100000000000000},
	};
	SharedCplexEnv env = LinearProgram::openCPLEX();
	std::vector<double> times;
	for (const auto& inst:instances) {
		clock_t start = std::clock();
		std::cout << "Solving instance " << inst.first << std::endl;
		std::ifstream in("special_tests/" + inst.first + ".tsp");
		if (!in) {
			std::cout << "Instance not found!" << std::endl;
			continue;
		}
		TSPInstance tsp(in);
		TSPSolution sol = tspsolvers::solveLP(tsp, nullptr, env, false);
		clock_t end = std::clock();
		double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
		if (sol.getCost() != inst.second) {
			std::cerr << "Failed to solve " << inst.first << ", found solution has cost " << sol.getCost()
					  << ", but optimal cost is " << inst.second << std::endl;
			break;
		} else {
			std::cout << "Found correct solution for " << inst.first << " in " << elapsed_secs << " seconds"
					  << std::endl;
		}
		times.push_back(elapsed_secs);
	}
	if (instances.size() == times.size()) {
		for (size_t i = 0; i < instances.size(); ++i) {
			std::cout << instances[i].first << ": " << times[i] << std::endl;
		}
	}
}
