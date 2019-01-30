#ifndef GRAPH_HPP
#define GRAPH_HPP
#include <istream>
#include <ostream>
#include <vector>
#include <cstdint>
#include <string>
#include "tour.hpp"


struct edge;


class Graph {
public:
	using node_id = unsigned;
	using edge_id = unsigned;
	using cost_t = unsigned;
	using Vec2 = std::array<double, 2>;
	using node = std::vector<edge_id>;

	explicit Graph(std::istream& input);

	explicit Graph(node_id nodeCount, const std::string& name);
	edge getEdge(edge_id id) const;
	edge_id addEdge(node_id endA, node_id endB, cost_t c = 0);
	Tour getReferenceTour();
	Tour greedyTSP() const;
	node_id getNodeCount() const;
	edge_id getEdgeCount() const;
private:
	void readNodes(std::istream& input, EdgeWeightType type);
	std::string name;
	std::vector<node> nodes;
	std::vector<edge> edges;
	std::vector<Vec2> nodeLocations;
	EdgeWeightType edgeType;
};

struct edge {
	Graph::node_id endA;
	Graph::node_id endB;
	Graph::cost_t cost;

	Graph::node_id getOtherNode(Graph::node_id known) const;
};

std::ostream& operator<<(std::ostream& out, const edge& e);

#endif
