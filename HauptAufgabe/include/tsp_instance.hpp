#ifndef TSP_INSTANCE_HPP
#define TSP_INSTANCE_HPP

#include <istream>
#include <vector>
#include <lemon/smart_graph.h>
#include <linear_program.hpp>

using cost_t = unsigned;
using node_id = int;
using edge_id = int;
using Graph = lemon::SmartGraph;

class TSPInstance {
public:
	//TODO are edge IDs consecutive? I can't find docs on that anywhere...
	explicit TSPInstance(std::istream& in);

	cost_t getDistance(node_id a, node_id b) const;

	node_id getSize() const;

	const Graph& getGraph();

	const Graph::EdgeMap <cost_t>& getGraphDistances();

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

	void setDistance(node_id a, node_id b, cost_t dist);

	Graph graph;
	Graph::EdgeMap <cost_t> graphDists;
	std::vector<std::vector<cost_t>> distances;
	std::string name;
};


#endif