#include <comb_heuristic.hpp>
#include <blossom_finder.hpp>

CombHeuristic::CombHeuristic(std::vector<ContractionRule *> rules) : rules(std::move(rules)), tolerance(1e-5) {}

std::vector<CombHeuristic::Comb> CombHeuristic::findViolatedCombs(const TspLpData& lpData, const TSPInstance& inst,
																  const std::vector<double>& sol) const {
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
	std::vector<Comb> ret;
	bool changed;
	do {
		changed = false;
		Graph::NodeMap<bool> used(g, false);
		for (ContractionRule *cr:rules) {
			changed |= cr->contractAll(g, used, capacity, contrMap);
		}
		if (changed) {
			BlossomFinder bf(g, capacity, tolerance, true);//TODO geht contractPaths hier?
			std::vector<BlossomFinder::Blossom> blossoms = bf.findViolatedBlossoms();
			for (const BlossomFinder::Blossom& b:blossoms) {
				if (b.isProperBlossom()) {
					ret.emplace_back(g, b, contrMap);
				}
			}
		}
	} while (changed && ret.size() < 60);
	return ret;
}

CombHeuristic::Comb::Comb(const Graph& g, const BlossomFinder::Blossom& b,
						  const tsp_util::ContractionMapTSP& contr) {
	for (Graph::Node n:b.handle) {
		handle.insert(handle.end(), contr[n].begin(), contr[n].end());
	}
	for (Graph::Edge e:b.teeth) {
		std::vector<city_id> tooth(contr[g.u(e)]);
		tooth.insert(tooth.end(), contr[g.v(e)].begin(), contr[g.v(e)].end());
		teeth.push_back(tooth);
	}
}

bool CombHeuristic::Comb::isBlossom() const {
	for (const std::vector<city_id>& tooth:teeth) {
		if (tooth.size() > 2) {
			return false;
		}
	}
	return true;
}
