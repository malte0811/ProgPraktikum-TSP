#include <contraction_rule.hpp>
#include <lemon/core.h>
#include <lemon/tolerance.h>
#include <utility>
#include <tsp_instance.hpp>
#include <tsp_utils.hpp>

ContractionRule::ContractionRule(ContractionRule::Validator rule)
		: validate(std::move(rule)) {}

bool ContractionRule::contractAll(Graph& g, lemon::GraphExtender<lemon::ListGraphBase>::NodeMap<bool>& used,
								  Graph::EdgeMap<double>& costs, tsp_util::ContractionMapTSP& contrMap) const {
	lemon::Tolerance<double> tolerance;
	bool changed = true;
	bool ret = false;
	do {
		//Alle noch nicht genutzten 1-Kanten und Knoten auflisten
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
		//Zu kontrahierende Mangen finden
		Contraction contr = validate(g, possibleNodes, possibleOneEdges, costs, used);
		if (!contr.empty()) {
			for (const std::vector<Graph::Node>& newNode:contr) {
				for (Graph::Node n:newNode) {
					used[n] = true;
				}
			}
			if (!contract(contr, g, costs, contrMap)) {
				changed = false;
			}
		} else {
			changed = false;
		}
		if (changed) {
			ret = true;
		}
	} while (changed);
	return ret;
}

bool ContractionRule::contract(const Contraction& contr, Graph& g, Graph::EdgeMap<double>& costs,
							   tsp_util::ContractionMapTSP& map) const {
	for (std::vector<Graph::Node> toContract:contr) {
		//Nur ein einzelner Knoten->nur enthalten, um als genutzt markiert zu werden
		if (toContract.size() < 2) {
			continue;
		}
		//Knoten, der nach der Kontraktion der Gesamtmenge entspricht
		Graph::Node remaining = toContract.back();
		toContract.pop_back();
		//Kosten der am verbleibenden Knoten anliegenden Kanten
		Graph::NodeMap<double> adjCosts(g, 0);
		for (Graph::OutArcIt it(g, remaining); it != lemon::INVALID; ++it) {
			adjCosts[g.target(it)] = costs[it];
		}
		for (Graph::Node toRemove:toContract) {
			for (Graph::OutArcIt it(g, toRemove); it != lemon::INVALID; ++it) {
				adjCosts[g.target(it)] += costs[it];
				if (adjCosts[g.target(it)] > 1) {
					return false;
				}
			}
		}
		std::vector<city_id>& contracted = map[remaining];
		for (Graph::Node toRemove:toContract) {
			contracted.insert(contracted.end(), map[toRemove].begin(), map[toRemove].end());
			g.erase(toRemove);
		}
		//Kosten der bereits existierenden Kanten erhöhen
		for (Graph::OutArcIt it(g, remaining); it != lemon::INVALID; ++it) {
			costs[it] = adjCosts[g.target(it)];
			adjCosts[g.target(it)] = 0;
		}
		//Neue Kanten hinzufügen (falls nötig)
		for (Graph::NodeIt it(g); it != lemon::INVALID; ++it) {
			if (adjCosts[it] > 0 && it != remaining) {
				Graph::Edge newE = g.addEdge(it, remaining);
				costs[newE] = adjCosts[it];
			}
		}
	}
	return true;
}
