#ifndef TSP_UTILS_HPP
#define TSP_UTILS_HPP

#include <linear_program.hpp>
#include <tsp_instance.hpp>

namespace tsp_util {
	std::vector<variable_id> createFractionalGraph(const TSPInstance& tsp, lemon::Tolerance<double> tolerance,
												   const std::vector<double>& solution, Graph& workGraph,
												   Graph::NodeMap <city_id>& workToOrig,
												   std::vector<Graph::Node>& origToWork,
												   Graph::EdgeMap <variable_id>& toVariable, Graph::EdgeMap<double>& c,
												   Graph::NodeMap<bool>& odd);
}


#endif