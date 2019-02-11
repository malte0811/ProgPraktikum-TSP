#ifndef TSP_INSTANCE_HPP
#define TSP_INSTANCE_HPP

#include <istream>
#include <vector>
#include <lemon/list_graph.h>
#include <linear_program.hpp>

using cost_t = unsigned;
using city_id = int;
using Graph = lemon::ListGraph;

class TSPInstance {
public:
	explicit TSPInstance(std::istream& in);

	inline cost_t getDistance(city_id a, city_id b) const;

	inline city_id getSize() const;

	inline const Graph& getGraph() const;

	inline const Graph::EdgeMap <cost_t>& getGraphDistances() const;

	inline variable_id getVariable(const Graph::Edge& e) const;

	inline variable_id getVariable(city_id a, city_id b) const;

	inline Graph::Edge getEdge(variable_id variable) const;

	inline city_id getCity(const Graph::Node& e) const;

	inline Graph::Node getNode(city_id city) const;

	inline variable_id getEdgeCount() const;

	void setupBasicLP(LinearProgram& lp) const;

	std::string getName() const;

	static constexpr city_id invalid_city = -1;

private:
	enum EdgeWeightType {
		euc_2d,
		ceil_2d,
		explicit_
	};
	enum EdgeFormat {
		full_matrix,
		lower_diag_row,
		upper_diag_row,
		upper_row
	};

	void readNodes(std::istream& input, EdgeWeightType type);

	void readEdges(std::istream& input, EdgeFormat type);

	void setDistance(city_id a, city_id b, cost_t dist);

	Graph graph;
	Graph::EdgeMap <cost_t> graphDists;
	Graph::EdgeMap <variable_id> edgeToVar;
	std::vector<Graph::Edge> varToEdge;
	Graph::NodeMap <city_id> nodeToCity;
	std::vector<Graph::Node> cityToNode;
	std::vector<std::vector<cost_t>> distances;
	std::string name;
};


city_id TSPInstance::getSize() const {
	return static_cast<city_id>(distances.size()+1);
}

const Graph& TSPInstance::getGraph() const {
	return graph;
}

const Graph::EdgeMap <cost_t>& TSPInstance::getGraphDistances() const {
	return graphDists;
}

variable_id TSPInstance::getVariable(const Graph::Edge& e) const {
	return edgeToVar[e];
}

variable_id TSPInstance::getVariable(city_id a, city_id b) const {
	if (b>a) {
		std::swap(a, b);
	}
	return (a*(a-1))/2+b;
}

Graph::Edge TSPInstance::getEdge(variable_id variable) const {
	return varToEdge[variable];
}

city_id TSPInstance::getCity(const Graph::Node& e) const {
	return nodeToCity[e];
}

Graph::Node TSPInstance::getNode(city_id city) const {
	return cityToNode[city];
}

variable_id TSPInstance::getEdgeCount() const {
	city_id size = getSize();
	return (size*(size-1))/2;
}

cost_t TSPInstance::getDistance(city_id a, city_id b) const {
	if (a>b) {
		return distances[a-1][b];
	} else if (b>a) {
		return distances[b-1][a];
	} else {
		return 0;
	}
}

#endif