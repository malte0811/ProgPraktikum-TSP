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
	int status;
	CPXENVptr env = CPXopenCPLEX(&status);
	if (status!=0) {
		throw std::runtime_error("Failed to open CPLEX environment: "+std::to_string(status));
	}
	LinearProgram lp(env, "assignment", LinearProgram::minimize);
	std::vector<double> min{0};
	std::vector<double> max{1};
	for (unsigned i = 0; i<=g.maxArcId(); ++i) {
		std::vector<double> obj{static_cast<double>(costs[lemon::SmartDigraph::arcFromId(i)])};
		lp.addVariables(obj, min, max);
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
			lp.addConstraint(<#initializer#>);
		} else if (nIt==end) {
			lp.addConstraint(<#initializer#>);
		} else {
			lp.addConstraint(<#initializer#>);
		}
	}
	BranchAndCut bb(lp, std::vector<CutGenerator*>(), 0);
	std::vector<long> sol = bb.solve();
	int totalCost = 0;
	for (int i = 0; i<sol.size(); ++i) {
		if (sol[i]>0) {
			auto arc = lemon::SmartDigraph::arcFromId(i);
			totalCost += costs[arc];
		}
	}
	CPXcloseCPLEX(&env);
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