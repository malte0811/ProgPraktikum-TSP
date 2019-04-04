#ifndef TSP_SOLVERS_HPP
#define TSP_SOLVERS_HPP

#include <tsp_instance.hpp>
#include <tsp_solution.hpp>
#include "tsp_lp_data.hpp"

namespace tspsolvers {
	TSPSolution solveGreedy(const TSPInstance& inst);

	TSPSolution solveLP(const TSPInstance& inst, const TSPSolution* initial, CPXENVptr& lpEnv, size_t maxOpenSize);

	void closeHamiltonPath(const TSPInstance& instance, std::vector<bool>& used, const std::vector<size_t>& degree,
						   const TspLpData& data);
}
#endif
