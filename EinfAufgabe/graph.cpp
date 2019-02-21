#include <iostream>
#include <ctime>
#include <istream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <stack>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include "graph.hpp"
#include "union_find.hpp"
#include "tour.hpp"

/*
Versucht, ein Objekt vom Typ T aus input zu lesen, und wirft
std::invalid_argument, falls dies nicht möglich war
*/
template<typename T>
T readOrThrow(std::istream& input) {
	T ret;
	input >> ret;
	if (!input) {
		throw std::invalid_argument("Could not read input!");
	}
	return ret;
}
Graph::Graph(std::istream& input) {
	std::string line;
	edgeType = euc_2d;
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
			std::cout << name << std::endl;
		} else if (keyword=="COMMENT") {
			//NOP
		} else if (keyword=="DIMENSION") {
			auto nodeCount = readOrThrow<node_id>(ss);
			nodes.resize(nodeCount);
			nodeLocations.resize(nodeCount);
		} else if (keyword=="EDGE_WEIGHT_TYPE") {
			auto type = readOrThrow<std::string>(ss);
			if (type=="EUC_2D") {
				edgeType = euc_2d;
			} else if (type=="CEIL_2D") {
				edgeType = ceil_2d;
			}
		} else if (keyword=="NODE_COORD_SECTION") {
			readNodes(input, edgeType);
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
	std::cout << "Finished reading" << std::endl;
}

void Graph::readNodes(std::istream& input, EdgeWeightType type) {
	std::vector<bool> set(getNodeCount(), false);
	node_id setCount = 0;
	for (node_id iteration = 0;iteration<getNodeCount();++iteration) {
		auto id = readOrThrow<node_id>(input)-1;
		assert(id<getNodeCount());
		assert(!set[id]);
		set[id] = true;
		++setCount;
		nodeLocations[id] = {readOrThrow<double>(input), readOrThrow<double>(input)};
	}
	assert(setCount==getNodeCount());
	for (node_id endA = 1;endA<getNodeCount();++endA) {
		for (node_id endB = 0;endB<endA;++endB) {
			Vec2 dist = nodeLocations[endA]-nodeLocations[endB];
			addEdge(endA, endB, dist.length(type));
		}
	}
}

Graph::Graph(node_id nodeCount, const std::string& name) : nodes(nodeCount), name(name), nodeLocations(nodeCount) {}

edge Graph::getEdge(edge_id id) const {
	return edges[id];
}

Graph::edge_id Graph::addEdge(node_id endA, node_id endB, cost_t c) {
	edge_id id = edges.size();
	edges.push_back({endA, endB, c});
	nodes[endA].push_back(id);
	nodes[endB].push_back(id);
	return id;
}

Graph::node_id Graph::getNodeCount() const {
	return nodes.size();
}

Graph::edge_id Graph::getEdgeCount() const {
	return edges.size();
}

const std::vector<Graph::edge_id>& Graph::getEdgesAt(Graph::node_id node) const {
	return nodes[node];
}

const std::string& Graph::getName() const {
	return name;
}

/*
Einfache Helfer-Methode, gibt den jeweils anderen Endknoten zu known zurück
*/
Graph::node_id edge::getOtherNode(Graph::node_id known) const {
	if (known==endA) {
		return endB;
	} else {
		return endA;
	}
}

std::ostream& operator<<(std::ostream& out, const edge& e) {
	out << "{" << e.endA << ", " << e.endB << "} with weight " << e.cost;
	return out;
}

Graph::cost_t Graph::Vec2::length(EdgeWeightType norm) const {
	double realLen = std::sqrt(x*x+y*y);
	switch (norm) {
		case euc_2d:
			return static_cast<cost_t>(std::round(realLen));
		case ceil_2d:
			return static_cast<cost_t>(std::ceil(realLen));
	}
	return 0;
}

Graph::Vec2 Graph::Vec2::operator-(const Graph::Vec2& other) const {
	return {x-other.x, y-other.y};
}
