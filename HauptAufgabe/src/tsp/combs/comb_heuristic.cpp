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
	} while (changed && ret.size() < 100);
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

size_t CombHeuristic::Comb::estimateNonzeroCount() const {
	size_t ret = (handle.size() * (handle.size() - 1)) / 2;
	for (const std::vector<city_id>& tooth:teeth) {
		ret += (tooth.size() * (tooth.size() - 1)) / 2;
	}
	return ret;
}

void CombHeuristic::Comb::simplify(const TspLpData& lpData, const std::vector<double>& solution) {
	const city_id invalid_city = TSPInstance::invalid_city;
	const city_id cityCount = lpData.getTSP().getCityCount();
	std::vector<bool> isHandle(cityCount);
	for (city_id i:handle) {
		isHandle[i] = true;
	}
	std::vector<bool> isPureHandle = isHandle;
	for (const auto& tooth:teeth) {
		for (city_id i:tooth) {
			isPureHandle[i] = false;
		}
	}
	lemon::Tolerance<double> tolerance;
	std::vector<std::pair<city_id, city_id>> neighbors(cityCount, {invalid_city, invalid_city});
	for (city_id i:handle) {
		if (!isPureHandle[i]) {
			continue;
		}
		std::vector<city_id> oneNeighbors;
		for (city_id neighbor = 0; neighbor < cityCount; ++neighbor) {
			if (!isPureHandle[neighbor]) {
				continue;
			}
			variable_id var = lpData.getVariable(i, neighbor);
			if (var != LinearProgram::invalid_variable) {
				double val = solution[var];
				if (!tolerance.different(1, val)) {
					//TODO bruache ich pure oder reicht handle?
					oneNeighbors.push_back(neighbor);
				}
			}
		}
		if (oneNeighbors.size() == 2) {
			neighbors[i] = {oneNeighbors[0], oneNeighbors[1]};
			//std::cout << "Found inner node " << i << " with neighbors " << oneNeighbors[0] << "," << oneNeighbors[1] << std::endl;
		} else if (oneNeighbors.size() == 1) {
			neighbors[i] = {oneNeighbors[0], invalid_city};
			//std::cout << "Found end node " << i << " with neighbor " << oneNeighbors[0] << std::endl;
		}
	}
	//std::cout << "Done" << std::endl;
	for (city_id i:handle) {
		std::pair<city_id, city_id> neighborsI = neighbors[i];
		if (neighborsI.first != invalid_city && neighborsI.second == invalid_city) {
			if (neighbors[neighborsI.first].second == invalid_city) {
				continue;
			}
			city_id last = i;
			city_id current = neighborsI.first;
			while (current != invalid_city) {
				if (last != i) {
					isHandle[last] = false;
				}
				city_id next;
				std::pair<city_id, city_id> neighborsCurr = neighbors[current];
				if (neighborsCurr.first != last) {
					next = neighborsCurr.first;
				} else {
					next = neighborsCurr.second;
				}
				last = current;
				current = next;
				assert(current != i);
			}
			city_id iOther = neighborsI.first;
			city_id lastOther = neighbors[last].first;
			assert(last != i);
			if (iOther == lastOther) {
				isHandle[iOther] = true;
			} else {
				assert(isHandle[i]);
				assert(isHandle[last]);
				assert(!isHandle[iOther]);
				assert(!isHandle[lastOther]);
				teeth.push_back({i, iOther});
				teeth.push_back({last, lastOther});
				neighbors[i] = neighbors[last] = {invalid_city, invalid_city};
			}
		}
	}
	for (size_t i = 0; i < handle.size();) {
		if (isHandle[handle[i]]) {
			++i;
		} else {
			handle[i] = handle.back();
			handle.pop_back();
		}
	}
	validate(cityCount);
}

void CombHeuristic::Comb::validate(city_id cityCount) const {
	std::vector<bool> isHandle(cityCount);
	for (city_id i:handle) {
		isHandle[i] = true;
	}
	std::vector<bool> inComb(cityCount);
	for (const auto& tooth:teeth) {
		bool foundHandle = false;
		bool foundNonHandle = false;
		for (city_id i:tooth) {
			assert(!inComb[i]);
			if (isHandle[i]) {
				foundHandle = true;
			} else {
				foundNonHandle = true;
			}
			inComb[i] = false;
		}
		assert(foundHandle);
		assert(foundNonHandle);
	}
}
