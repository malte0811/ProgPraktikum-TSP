#ifndef TSP_SOLVERS_HPP
#define TSP_SOLVERS_HPP

#include <tsp_instance.hpp>
#include <tsp_solution.hpp>
#include "tsp_lp_data.hpp"

namespace tspsolvers {
	namespace cutgens {
		extern const char *const connected;
		extern const char *const subtour;
		extern const char *const simpleCombs;
		extern const char *const twoMatching;
		extern const char *const generalCombs;
		extern const char *const defaultGens;
	}

	TSPSolution solveGreedy(const TSPInstance& inst);

	TSPSolution solveLP(const TSPInstance& inst, const TSPSolution *initial, CPXENVptr& lpEnv, size_t maxOpenSize,
						const std::vector<std::string>& cutGenerators =
								{cutgens::connected, cutgens::subtour, cutgens::simpleCombs,
								 cutgens::twoMatching, cutgens::generalCombs});

	void closeHamiltonPath(const lemon::FullGraph& g, lemon::FullGraph::EdgeMap<bool>& used,
						   const lemon::FullGraph::NodeMap<size_t>& degree);
}
#endif
