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
#include <tsp_solution.hpp>
#include <union_find.hpp>
#include <tsp_solvers.hpp>
#include <subtour_cut_gen.hpp>
#include <two_matching_cut_gen.hpp>
#include <branch_and_cut.hpp>

namespace tspsolvers {
	TSPSolution solveGreedy(const TSPInstance& inst) {
		struct sorting_data {
			Graph::Edge e;
			cost_t cost;
		};
		const Graph::EdgeMap <cost_t>& distances = inst.getGraphDistances();
		std::vector<sorting_data> sortedEdges;
		for (Graph::EdgeIt it(inst.getGraph()); it!=lemon::INVALID; ++it) {
			sortedEdges.push_back({it, distances[it]});
		}
		std::sort(sortedEdges.begin(), sortedEdges.end(),
				  [](const sorting_data& edgeA, const sorting_data& edgeB) {
					  return edgeA.cost<edgeB.cost;
				  }
		);
		UnionFind connectedComps(inst.getSize());
		std::vector<bool> used(inst.getEdgeCount());
		Graph::NodeMap <size_t> degree(inst.getGraph());
		unsigned addedEdges = 0;
		cost_t totalCost = 0;
		for (const sorting_data& data:sortedEdges) {
			const Graph::Edge& e = data.e;
			Graph::Node u = inst.getGraph().u(e);
			size_t& edgesAtU = degree[u];
			if (edgesAtU>=2) {
				continue;
			}
			Graph::Node v = inst.getGraph().v(e);
			size_t& edgesAtV = degree[v];
			if (edgesAtV>=2) {
				continue;
			}
			city_id rootU = connectedComps.find(inst.getCity(u));
			city_id rootV = connectedComps.find(inst.getCity(v));
			//Unterschiedliche Zusammenhangskomponenten->hinzufügen
			if (rootU!=rootV) {
				++addedEdges;
				++edgesAtU;
				++edgesAtV;
				used[inst.getVariable(e)] = true;
				if (addedEdges==inst.getSize()-1) {
					break;
				}
				connectedComps.mergeRoots(rootU, rootV);
			}
		}
		city_id firstEnd = TSPInstance::invalid_city;
		for (Graph::NodeIt it(inst.getGraph()); it!=lemon::INVALID; ++it) {
			if (degree[it]<2) {
				if (firstEnd==TSPInstance::invalid_city) {
					firstEnd = inst.getCity(it);
				} else {
					city_id otherEnd = inst.getCity(it);
					used[inst.getVariable(firstEnd, otherEnd)] = true;
					break;
				}
			}
		}
		return TSPSolution(inst, used);
	}

	TSPSolution solveLP(const TSPInstance& inst, const TSPSolution* initial) {
		LinearProgram lp(inst.getName(), LinearProgram::minimize);
		inst.setupBasicLP(lp);
		SubtourCutGen subtours(inst);
		TwoMatchingCutGen matchings(inst);
		BranchAndCut bac(lp, {&subtours, &matchings});
		if (initial!=nullptr) {
			std::vector<long> asVars(inst.getEdgeCount());
			for (variable_id i = 0; i<inst.getEdgeCount(); ++i) {
				asVars[i] = initial->getVariable(i);
			}
			bac.setUpperBound(asVars, initial->getCost());
		}
		std::vector<long> tour = bac.solve();
		std::vector<bool> asBools(inst.getEdgeCount());
		for (variable_id i = 0; i<inst.getEdgeCount(); ++i) {
			asBools[i] = tour[i];
		}
		return TSPSolution(inst, asBools);
	}
}