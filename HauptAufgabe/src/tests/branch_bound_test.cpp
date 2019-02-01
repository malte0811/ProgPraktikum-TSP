#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <linear_program.hpp>
#include <branch_and_cut.hpp>
#include <fstream>
#include <lemon/smart_graph.h>

int shortestPath(const std::string& name) {
	lemon::SmartDigraph g;
	lemon::SmartDigraph::ArcMap<int> costs(g);
	lemon::SmartDigraph::Node start = lemon::SmartDigraph::nodeFromId(0), end = lemon::SmartDigraph::nodeFromId(1);
	std::ifstream in(name);
	if (!in) {
		return -13;
	}
	unsigned nodeCount;
	in >> nodeCount >> std::ws;
	for (unsigned i = 0; i<nodeCount; ++i) {
		g.addNode();
	}
	std::string line;
	while (std::getline(in, line)) {
		std::stringstream ss(line);
		unsigned edgeA, edgeB;
		int cost;
		ss >> edgeA >> edgeB >> cost;
		if (!ss) break;
		auto e = g.addArc(lemon::SmartDigraph::nodeFromId(edgeA), lemon::SmartDigraph::nodeFromId(edgeB));
		costs[e] = cost;
	}
	LinearProgram lp("assignment", LinearProgram::minimize);
	for (unsigned i = 0; i<=g.maxArcId(); ++i) {
		lp.addVariable(costs[lemon::SmartDigraph::arcFromId(i)], 0, 1);
	}
	for (lemon::SmartDigraph::NodeIt nIt(g); nIt!=lemon::INVALID; ++nIt) {
		std::vector<int> indices;
		for (lemon::SmartDigraph::OutArcIt it(g, nIt); it!=lemon::INVALID; ++it) {
			indices.push_back(lemon::SmartDigraph::id(it));
		}
		std::vector<double> values(indices.size(), 1);
		for (lemon::SmartDigraph::InArcIt it(g, nIt); it!=lemon::INVALID; ++it) {
			indices.push_back(lemon::SmartDigraph::id(it));
		}
		values.resize(indices.size(), -1);
		if (nIt==start) {
			lp.addConstraint(indices, values, 1, LinearProgram::equal);
		} else if (nIt==end) {
			lp.addConstraint(indices, values, -1, LinearProgram::equal);
		} else {
			lp.addConstraint(indices, values, 0, LinearProgram::equal);
		}
	}
	BranchAndCut bb(lp, std::vector<CutGenerator*>());
	std::vector<long> sol = bb.solve();
	int totalCost = 0;
	for (int i = 0; i<sol.size(); ++i) {
		if (sol[i]>0) {
			auto arc = lemon::SmartDigraph::arcFromId(i);
			totalCost += costs[arc];
		}
	}
	return totalCost;
}

int main() {
	const unsigned len = 8;
	int expected[len] = {810, -1341, 1623, 1108, -600, 16, 580, 938};
	for (unsigned i = 0; i<len; ++i) {
		std::cout << shortestPath("data/conservative"+std::to_string(i+1)) << ", expected value: " << expected[i]
				  << std::endl;
	}
}