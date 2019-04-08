#ifndef TSP_UTILS_HPP
#define TSP_UTILS_HPP

#include <linear_program.hpp>
#include <tsp_instance.hpp>
#include "tsp_lp_data.hpp"

namespace tsp_util {
	using ContractionMapTSP = Graph::NodeMap<std::vector<city_id>>;
	using ContractionMapGraph = Graph::NodeMap<std::vector<Graph::Node>>;
	std::vector<variable_id> createFractionalGraph(const TSPInstance& tsp, const TspLpData& lpData,
												   lemon::Tolerance<double> tolerance,
												   const std::vector<double>& solution, Graph& workGraph,
												   Graph::NodeMap <city_id>& workToOrig,
												   std::vector<Graph::Node>& origToWork,
												   Graph::EdgeMap <variable_id>& toVariable,
												   Graph::EdgeMap<double>& c, Graph::NodeMap<bool>& odd);

	void addSupportGraphEdges(const TSPInstance& tsp, const TspLpData& lpData, lemon::Tolerance<double> tolerance,
							  const std::vector<double>& solution, Graph& workGraph,
							  std::vector<Graph::Node>& origToWork,
							  Graph::NodeMap <city_id>& workToOrig, Graph::EdgeMap<double>& c);
}


#endif