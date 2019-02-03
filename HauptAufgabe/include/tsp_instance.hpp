#ifndef TSP_INSTANCE_HPP
#define TSP_INSTANCE_HPP

#include <istream>
#include <vector>
#include <lemon/smart_graph.h>
#include <linear_program.hpp>

using cost_t = unsigned;
using city_id = int;
using Graph = lemon::SmartGraph;

class TSPInstance {
public:
	explicit TSPInstance(std::istream& in);

	cost_t getDistance(city_id a, city_id b) const;

	city_id getSize() const;

	const Graph& getGraph() const;

	const Graph::EdgeMap <cost_t>& getGraphDistances();

	variable_id getVariable(const Graph::Edge& e) const;

	Graph::Edge getEdge(variable_id variable) const;

	city_id getCity(const Graph::Node& e) const;

	Graph::Node getNode(city_id city) const;

	void setupBasicLP(LinearProgram& lp);

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


#endif