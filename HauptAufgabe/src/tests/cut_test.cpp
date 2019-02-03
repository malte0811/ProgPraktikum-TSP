#include <lemon/smart_graph.h>
#include <cut_generator.hpp>
#include <fstream>
#include <branch_and_cut.hpp>
#include <cmath>

class PathCutGenerator : public CutGenerator {
public:
	PathCutGenerator(const lemon::SmartDigraph& graph, lemon::SmartDigraph::Node start, lemon::SmartDigraph::Node end)
			: graph(graph), start(start), end(end) {}

	bool validate(LinearProgram& lp, const std::vector<double>& solution) override {
		for (lemon::SmartDigraph::NodeIt it(graph); it!=lemon::INVALID; ++it) {
			double sum = 0;
			std::vector<int> edges;
			for (lemon::SmartDigraph::OutArcIt eIt(graph, it); eIt!=lemon::INVALID; ++eIt) {
				sum += solution[lemon::SmartDigraph::id(eIt)];
				edges.push_back(lemon::SmartDigraph::id(eIt));
			}
			const unsigned outDegree = edges.size();
			std::vector<double> coeffs(outDegree, 1);
			for (lemon::SmartDigraph::InArcIt eIt(graph, it); eIt!=lemon::INVALID; ++eIt) {
				sum -= solution[lemon::SmartDigraph::id(eIt)];
				edges.push_back(lemon::SmartDigraph::id(eIt));
			}
			coeffs.resize(edges.size(), -1);
			//TODO rounding errors?
			if (it==start && std::fabs(sum-1)>eps) {
				lp.addConstraint(edges, coeffs, 1, LinearProgram::equal);
				std::cout << "Adding start cut" << std::endl;
				return false;
			} else if (it==end && std::fabs(sum+1)>eps) {
				lp.addConstraint(edges, coeffs, -1, LinearProgram::equal);
				std::cout << "Adding end cut" << std::endl;
				return false;
			} else if (it!=start && it!=end && std::fabs(sum)>eps) {
				lp.addConstraint(edges, coeffs, 0, LinearProgram::equal);
				int id = lemon::SmartDigraph::id(it);
				if (id%100==0) std::cout << "Adding cut for " << id << std::endl;
				return false;
			}
		}
		return true;//TODO
	}

private:
	static constexpr double eps = 0.001;
	const lemon::SmartDigraph& graph;
	lemon::SmartDigraph::Node start;
	lemon::SmartDigraph::Node end;
};

int shortestPath(const std::string& name) {
	lemon::SmartDigraph g;
	lemon::SmartDigraph::ArcMap<int> costs(g);
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
	//TODO is this necessary? Why?

	std::vector<int> indices;
	for (lemon::SmartDigraph::OutArcIt it(g, lemon::SmartDigraph::nodeFromId(0)); it!=lemon::INVALID; ++it) {
		indices.push_back(lemon::SmartDigraph::id(it));
	}
	std::vector<double> values(indices.size(), 1);
	for (lemon::SmartDigraph::InArcIt it(g, lemon::SmartDigraph::nodeFromId(0)); it!=lemon::INVALID; ++it) {
		indices.push_back(lemon::SmartDigraph::id(it));
	}
	values.resize(indices.size(), -1);
	lp.addConstraint(indices, values, 1, LinearProgram::equal);
	PathCutGenerator gen(g, lemon::SmartDigraph::nodeFromId(0), lemon::SmartDigraph::nodeFromId(1));
	BranchAndCut bb(lp, {&gen});
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

//Sehr langsam, gibt aber korrekte Ergebnisse. Kürzeste Wege sind halt ein dummer Testfall.
int main() {
	const unsigned len = 8;
	int expected[len] = {810, -1341, 1623, 1108, -600, 16, 580, 938};
	for (unsigned i = 0; i<len; ++i) {
		std::cout << shortestPath("data/conservative"+std::to_string(i+1)) << ", expected value: " << expected[i]
				  << std::endl;
	}
}