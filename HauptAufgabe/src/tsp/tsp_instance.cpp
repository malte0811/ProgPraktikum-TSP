#include <tsp_instance.hpp>
#include <sstream>
#include <cmath>
#include <array>
#include <cassert>

template<typename T>
T readOrThrow(std::istream& input) {
	T ret;
	input >> ret;
	if (!input) {
		throw std::invalid_argument("Could not read input!");
	}
	return ret;
}

TSPInstance::TSPInstance(std::istream& input) : graphDists(graph) {
	std::string line;
	EdgeWeightType edgeType = euc_2d;
	EdgeFormat edgeFormat = full_matrix;
	while (std::getline(input, line) && !line.empty()) {
		std::stringstream ss(line);
		std::string keyword;
		ss >> keyword >> std::ws;
		if (ss.peek()==':') {
			ss.ignore();
		} else if (keyword.back()==':') {
			keyword = keyword.substr(0, keyword.size()-1);
		}
		if (keyword=="NAME") {
			ss >> name;
		} else if (keyword=="COMMENT") {
			//NOP
		} else if (keyword=="DIMENSION") {
			auto nodeCount = readOrThrow<node_id>(ss);
			distances.resize(nodeCount-1);
			for (node_id i = 0; i<nodeCount-1; ++i) {
				distances[i].resize(i+1);
			}
		} else if (keyword=="EDGE_WEIGHT_TYPE") {
			auto type = readOrThrow<std::string>(ss);
			if (type=="EUC_2D") {
				edgeType = euc_2d;
			} else if (type=="CEIL_2D") {
				edgeType = ceil_2d;
			} else if (type=="EXPLICIT") {
				edgeType = explicit_;
			} else {
				throw std::runtime_error("Unknown edge weight type: "+type);
			}
		} else if (keyword=="EDGE_WEIGHT_FORMAT") {
			auto type = readOrThrow<std::string>(ss);
			if (type=="FULL_MATRIX") {
				edgeFormat = full_matrix;
			} else if (type=="LOWER_DIAG_ROW") {
				edgeFormat = lower_diag_row;
			} else if (type=="UPPER_DIAG_ROW") {
				edgeFormat = upper_diag_row;
			} else if (type=="UPPER_ROW") {
				edgeFormat = upper_row;
			} else {
				throw std::runtime_error("Unknown edge format: "+type);
			}
		} else if (keyword=="NODE_COORD_SECTION") {
			readNodes(input, edgeType);
		} else if (keyword=="EDGE_WEIGHT_SECTION") {
			readEdges(input, edgeFormat);
		} else if (keyword=="TYPE") {
			auto type = readOrThrow<std::string>(ss);
			if (type!="TSP") {
				throw std::invalid_argument("Input is not a symmetric TSP instance!");
			}
		} else if (keyword=="EOF") {
			break;
		} else {
			throw std::invalid_argument("Invalid keyword in input: \""+keyword+"\"");
		}
	}
	node_id nodeCount = getSize();
	for (node_id id = 0; id<nodeCount; ++id) {
		Graph::Node nId = graph.addNode();
		for (node_id other = 0; other<id; ++other) {
			Graph::Edge e = graph.addEdge(nId, Graph::nodeFromId(other));
			graphDists[e] = getDistance(id, other);
		}
	}
}

void TSPInstance::readNodes(std::istream& input, EdgeWeightType type) {
	node_id nodeCount = getSize();
	std::vector<std::array<double, 2>> nodeLocations(nodeCount);
	for (node_id iteration = 0; iteration<nodeCount; ++iteration) {
		auto id = readOrThrow<node_id>(input)-1;
		assert(id<nodeCount);
		nodeLocations[id][0] = readOrThrow<double>(input);
		nodeLocations[id][1] = readOrThrow<double>(input);
	}
	//TODO check that all locations have been set
	for (node_id higherId = 1; higherId<nodeCount; ++higherId) {
		for (node_id lowerId = 0; lowerId<higherId; ++lowerId) {
			cost_t distance = 0;
			double distX = nodeLocations[lowerId][0]-nodeLocations[higherId][0];
			double distY = nodeLocations[lowerId][1]-nodeLocations[higherId][1];
			switch (type) {
				case euc_2d:
					distance = std::round(std::sqrt(distX*distX+distY*distY));
					break;
				case ceil_2d:
					distance = std::ceil(std::sqrt(distX*distX+distY*distY));
					break;
				case explicit_:
					throw std::runtime_error("Weight type is EXPLICIT, but a NODE_COORD_SECTION exists!");
			}
			distances[higherId-1][lowerId] = distance;
		}
	}
}

void TSPInstance::readEdges(std::istream& input, TSPInstance::EdgeFormat type) {
	const node_id nodeCount = getSize();
	node_id rowCount = nodeCount;
	if (type==upper_row) {
		rowCount--;
	}
	for (node_id row = 0; row<rowCount; ++row) {
		node_id minCol = 0;
		node_id colCount = 0;
		switch (type) {
			case full_matrix:
				colCount = nodeCount;
				break;
			case lower_diag_row:
				colCount = row+1;
				break;
			case upper_diag_row:
				minCol = row;
				colCount = nodeCount-row;
				break;
			case upper_row:
				minCol = row+1;
				colCount = nodeCount-row+1;
				break;
		}
		for (node_id col = minCol; col<minCol+colCount; ++col) {
			setDistance(row, col, readOrThrow<cost_t>(input));
		}
	}
}

cost_t TSPInstance::getDistance(node_id a, node_id b) const {
	if (a>b) {
		return distances[a-1][b];
	} else if (b>a) {
		return distances[b-1][a];
	} else {
		return 0;
	}
}

void TSPInstance::setDistance(node_id a, node_id b, cost_t dist) {
	if (a>b) {
		distances[a-1][b] = dist;
	} else if (b>a) {
		distances[b-1][a] = dist;
	}
}

node_id TSPInstance::getSize() const {
	return distances.size()+1;
}

const Graph& TSPInstance::getGraph() {
	return graph;
}

const Graph::EdgeMap <cost_t>& TSPInstance::getGraphDistances() {
	return graphDists;
}

void TSPInstance::setupBasicLP(LinearProgram& lp) {
	for (edge_id i = 0; i<=graph.maxEdgeId(); ++i) {
		Graph::Edge e = Graph::edgeFromId(i);
		lp.addVariable(graphDists[e], 0, 1);
	}
	for (Graph::NodeIt nIt(graph); nIt!=lemon::INVALID; ++nIt) {
		std::vector<edge_id> adjancent;
		for (Graph::IncEdgeIt eIt(graph, nIt); eIt!=lemon::INVALID; ++eIt) {
			adjancent.push_back(Graph::id(eIt));
		}
		lp.addConstraint(adjancent, std::vector<double>(adjancent.size(), 1), 2, LinearProgram::equal);
	}
}
