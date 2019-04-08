#include "comb_heuristic.hpp"

CombHeuristic::CombHeuristic(std::vector<ContractionRule *> rules) : rules(std::move(rules)), tolerance(1e-5) {}

std::vector<CombHeuristic::Comb> CombHeuristic::findViolatedCombs(const TspLpData& lpData, const TSPInstance& inst,
																  const std::vector<double>& sol) {
	Graph g;
	std::vector<Graph::Node> origToWork;
	Graph::EdgeMap<double> capacity(g);
	tsp_util::ContractionMapTSP contrMap(g);
	{
		Graph::NodeMap <city_id> workToOrig(g);
		tsp_util::addSupportGraphEdges(inst, lpData, tolerance, sol, g, origToWork, workToOrig, capacity);
		for (Graph::NodeIt it(g); it != lemon::INVALID; ++it) {
			contrMap[it] = {workToOrig[it]};
		}
	}
	bool changed;
	do {
		changed = false;
		Graph::NodeMap<bool> used(g, false);
		for (ContractionRule *cr:rules) {
			changed |= cr->contractAll(g, used, capacity, contrMap);
		}

	} while (changed);
}
