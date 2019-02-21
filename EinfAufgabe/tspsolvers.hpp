#ifndef GREEDY_TSP_HPP
#define GREEDY_TSP_HPP

#include "tour.hpp"
#include "graph.hpp"

namespace tspsolvers {
	Tour greedyTSP(const Graph& g);

	Tour getReferenceTour(const Graph& g);

	std::vector<Graph::node_id> getOrderFromEdges(const std::vector<std::vector<edge>>& edgesAtNode);

	void closeHamiltonPath(const Graph& g, std::vector<std::vector<edge>>& edgesAtNode);
}
#endif