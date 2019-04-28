#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <istream>
#include <ostream>
#include <vector>
#include <cstdint>
#include <string>
#include <array>

struct edge;

enum EdgeWeightType {
	euc_2d,
	ceil_2d
};

class Graph {
public:
	using node_id = unsigned;
	using edge_id = unsigned;
	using cost_t = unsigned;

	struct Vec2 {
		double x, y;

		cost_t length(EdgeWeightType norm) const;

		Vec2 operator-(const Vec2& other) const;
	};

	using node = std::vector<edge_id>;

	explicit Graph(std::istream& input);

	edge getEdge(edge_id id) const;

	edge_id addEdge(node_id endA, node_id endB, cost_t c = 0);

	node_id getNodeCount() const;

	edge_id getEdgeCount() const;

	const std::vector<edge_id>& getEdgesAt(node_id node) const;

	const std::string& getName() const;

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
