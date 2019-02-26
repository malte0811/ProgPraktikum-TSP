#ifndef TSP_SOLVERS_HPP
#define TSP_SOLVERS_HPP

#include <tsp_instance.hpp>
#include <tsp_solution.hpp>

namespace tspsolvers {
	TSPSolution solveGreedy(const TSPInstance& inst);

	TSPSolution solveLP(const TSPInstance& inst, const TSPSolution* initial, CPXENVptr& lpEnv, size_t maxOpenSize);

	void closeHamiltonPath(const TSPInstance& instance, std::vector<bool>& used, const Graph::NodeMap <size_t>& degree);
}
#endif
