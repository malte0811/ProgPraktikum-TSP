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
	const std::vector<std::string> instances{
			"data_and_type",
			"data_and_type_2",
			"duplicate_data",
			"duplicate_keys",
			"empties",
			"empty_nodes",
			"ewsection_function",
			"few_nodes",
			"inv_edge_format",
			"inv_edge_type",
			"many_nodes",
			"many_nodes_2",
			"no_distance",
			"no_name",
			"nodec_before_dim",
			"non_tsp",
	};
	for (const std::string& inst:instances) {
		std::ifstream in("format_tests/" + inst + ".tsp");
		if (!in) {
			std::cout << "Instance not found!" << std::endl;
			continue;
		}
		try {
			TSPInstance tsp(in);
			std::cout << inst << " did not throw an error!" << std::endl;
		} catch (const std::runtime_error& err) {
			std::cout << inst << ": " << err.what() << std::endl;
		}
	}
}
