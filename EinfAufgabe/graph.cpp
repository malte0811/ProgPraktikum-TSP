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
	for (node_id iteration = 0;iteration<getNodeCount();++iteration) {
		auto id = readOrThrow<node_id>(input)-1;
		assert(id<getNodeCount());
		nodeLocations[id][0] = readOrThrow<double>(input);
		nodeLocations[id][1] = readOrThrow<double>(input);
	}
	//TODO check that all locations have been set
	for (node_id endA = 1;endA<getNodeCount();++endA) {
		for (node_id endB = 0;endB<endA;++endB) {
			cost_t distance = 0;
			double distX = nodeLocations[endB][0]-nodeLocations[endA][0];
			double distY = nodeLocations[endB][1]-nodeLocations[endA][1];
			switch (type) {
				case euc_2d:
					//TODO
					distance = std::round(std::sqrt(distX*distX+distY*distY));
					break;
				case ceil_2d:
					distance = std::ceil(std::sqrt(distX*distX+distY*distY));
					break;
			}
			addEdge(endA, endB, distance);
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

Tour Graph::greedyTSP() const {
	/*
	sortedEdges enthält sorting_data-Objekte, so dass die zu den ID's gehörenden
	Kanten aufsteigend nach ihren Kosten sortiert sind. Es ist schneller, ein
	struct zu verwenden, als nur die IDs zu speichern und beim Sortieren auf edges
	zuzugreifen.
	*/
	typedef struct {
		edge_id id;
		cost_t cost;
	} sorting_data;
	std::vector<sorting_data> sortedEdges(edges.size());
	for (edge_id i = 0;i<edges.size();++i) {
		sortedEdges[i] = {i, edges[i].cost};
	}
	//Sortieren nach Kosten der entsprechenden Kanten
	std::sort(sortedEdges.begin(), sortedEdges.end(),
		[this](const sorting_data& edgeA, const sorting_data& edgeB) {
			return edgeA.cost<edgeB.cost;
		}
	);
	//Union-Find-Struktur, in der die Zusammenhangskomponenten gespeichert werden
	UnionFind connectedComps(getNodeCount());
	//TODO nicer structure?
	std::vector<std::vector<edge_id>> edgesAtNode(getNodeCount());
	unsigned addedEdges = 0;
	cost_t totalCost = 0;
	for (const sorting_data& data:sortedEdges) {
		const edge& e = edges[data.id];
		std::vector<edge_id>& edgesAtA = edgesAtNode[e.endA];
		if (edgesAtA.size()>=2) {
			continue;
		}
		std::vector<edge_id>& edgesAtB = edgesAtNode[e.endB];
		if (edgesAtB.size()>=2) {
			continue;
		}
		node_id rootA = connectedComps.find(e.endA);
		node_id rootB = connectedComps.find(e.endB);
		//Unterschiedliche Zusammenhangskomponenten->hinzufügen
		if (rootA!=rootB || addedEdges==nodes.size()-1) {
			++addedEdges;
			totalCost += e.cost;
			edgesAtA.push_back(data.id);
			edgesAtB.push_back(data.id);
			if (addedEdges==nodes.size()) {
				break;//MST ist vollständig
			}
			//mergeRoots, um find zu überspringen
			connectedComps.mergeRoots(rootA, rootB);
		}
	}
	std::vector<node_id> order;
	order.reserve(getNodeCount());
	node_id prevNode = 0;
	node_id currentNode = 0;
	do {
		for (edge_id eid:edgesAtNode[currentNode]) {
			const edge& e = getEdge(eid);
			node_id otherNode = e.getOtherNode(currentNode);
			if (otherNode!=prevNode) {
				order.push_back(currentNode);
				prevNode = currentNode;
				currentNode = otherNode;
				break;
			}
		}
	} while (currentNode!=0);
	return Tour(nodeLocations, order, edgeType, totalCost);
}

Tour Graph::getReferenceTour() {
	std::vector<node_id> order;
	order.reserve(getNodeCount());
	cost_t totalCost = 0;
	for (node_id curr = 0;curr<getNodeCount();++curr) {
		order.push_back(curr);
		for (edge_id eid:nodes[curr]) {
			const edge& e = getEdge(eid);
			if (e.getOtherNode(curr)==(curr+1)%getNodeCount()) {
				totalCost += e.cost;
				break;
			}
		}
	}
	return Tour(nodeLocations, order, edgeType, totalCost);
}

Graph::node_id Graph::getNodeCount() const {
	return nodes.size();
}

Graph::edge_id Graph::getEdgeCount() const {
	return edges.size();
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
