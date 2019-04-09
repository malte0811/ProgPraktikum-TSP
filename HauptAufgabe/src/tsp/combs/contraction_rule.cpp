#include <contraction_rule.hpp>

ContractionRule::ContractionRule(ContractionRule::Validator rule)
		: validate(std::move(rule)) {}

bool ContractionRule::contractAll(Graph& g, lemon::GraphExtender<lemon::ListGraphBase>::NodeMap<bool>& used,
								  Graph::EdgeMap<double>& costs, tsp_util::ContractionMapTSP& contrMap) const {
	lemon::Tolerance<double> tolerance;
	bool changed = true;
	bool ret = false;
	do {
		std::vector<Graph::Node> possibleNodes;
		std::vector<Graph::Edge> possibleOneEdges;
		for (Graph::NodeIt it(g); it != lemon::INVALID; ++it) {
			if (!used[it]) {
				possibleNodes.push_back(it);
			}
		}
		for (Graph::EdgeIt it(g); it != lemon::INVALID; ++it) {
			if (!tolerance.different(1, costs[it]) && !used[g.u(it)] && !used[g.v(it)]) {
				possibleOneEdges.push_back(it);
			}
		}
		Contraction contr = validate(g, possibleNodes, possibleOneEdges, costs, used);
		if (!contr.empty()) {
			for (const std::vector<Graph::Node>& newNode:contr) {
				for (Graph::Node n:newNode) {
					used[n] = true;
				}
			}
			contract(contr, g, costs, contrMap);
		} else {
			changed = false;
		}
		if (changed) {
			ret = true;
		}
	} while (changed);
	return ret;
}

void ContractionRule::contract(const Contraction& contr, Graph& g, Graph::EdgeMap<double>& costs,
							   tsp_util::ContractionMapTSP& map) const {
	Graph::NodeMap<bool> used(g, false);
	for (std::vector<Graph::Node> toContract:contr) {
		for (auto n:toContract) {
			assert(!used[n]);
			assert(g.valid(n));
			used[n] = true;
		}
		if (toContract.size() < 2) {
			continue;
		}
		Graph::Node remaining = toContract.back();
		toContract.pop_back();
		Graph::NodeMap<double> adjCosts(g, 0);
		for (Graph::OutArcIt it(g, remaining); it != lemon::INVALID; ++it) {
			adjCosts[g.target(it)] = costs[it];
		}
		std::vector<city_id>& contracted = map[remaining];
		for (Graph::Node toRemove:toContract) {
			contracted.insert(contracted.end(), map[toRemove].begin(), map[toRemove].end());
			for (Graph::OutArcIt it(g, toRemove); it != lemon::INVALID; ++it) {
				adjCosts[g.target(it)] += costs[it];
			}
			g.erase(toRemove);
		}
		for (Graph::OutArcIt it(g, remaining); it != lemon::INVALID; ++it) {
			costs[it] = adjCosts[g.target(it)];
			adjCosts[g.target(it)] = 0;
		}
		for (Graph::NodeIt it(g); it != lemon::INVALID; ++it) {
			if (adjCosts[it] > 0 && it != remaining) {
				Graph::Edge newE = g.addEdge(it, remaining);
				costs[newE] = adjCosts[it];
			}
		}
	}
}
