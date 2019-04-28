#include <comb_heuristic.hpp>
#include <blossom_finder.hpp>

CombHeuristic::CombHeuristic(std::vector<ContractionRule *> rules) : rules(std::move(rules)), tolerance(1e-5) {}

std::vector<CombHeuristic::Comb> CombHeuristic::findViolatedCombs(const TspLpData& lpData, const TSPInstance& inst,
																  const std::vector<double>& sol) const {
	Graph g;
	std::vector<Graph::Node> origToWork;
	Graph::EdgeMap<double> capacity(g);
	tsp_util::ContractionMapTSP contrMap(g);
	//Graph und Kontraktionsdaten initialisieren
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
		//Graph nach den vorgegebenen Mustern kontrahieren, ohne Knoten doppelt zu nutzen
		Graph::NodeMap<bool> used(g, false);
		for (ContractionRule *cr:rules) {
			changed |= cr->contractAll(g, used, capacity, contrMap);
		}
		if (changed) {
			//Verletzte Blüten im kontrahierten Graphen finden
			BlossomFinder bf(g, capacity, tolerance, true);
			std::vector<BlossomFinder::Blossom> blossoms = bf.findViolatedBlossoms();
			for (const BlossomFinder::Blossom& b:blossoms) {
				if (b.isProperBlossom()) {
					ret.emplace_back(g, b, contrMap);
				}
			}
		}
		//Abbrechen, wenn keine Mengen mehr kontrahiert wurden oder mehr als 100 Combs gefunden sind
	} while (changed && ret.size() < 100);
	return ret;
}

//Wandelt eine Blüte im kontrahierten Graphen in einen Kamm im TSP-Graphen um
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

//Wahr, falls der Kamm eine Blüte ist, d.h. alle Zinken Größe 2 haben
bool CombHeuristic::Comb::isBlossom() const {
	for (const std::vector<city_id>& tooth:teeth) {
		if (tooth.size() > 2) {
			return false;
		}
	}
	return true;
}

//Vereinfacht den Kamm, z.B. durch das Entfernen langer Ketten von 1-Kanten aus dem Griff
void CombHeuristic::Comb::simplify(const TspLpData& lpData, const std::vector<double>& solution) {
	const city_id cityCount = lpData.getTSP().getCityCount();
	//Ist die Stadt im Griff enthalten?
	std::vector<bool> isHandle(cityCount);
	for (city_id i:handle) {
		isHandle[i] = true;
	}
	//Ist die Stadt im Griff, aber in keinem Zahn enthalten?
	std::vector<bool> isPureHandle = isHandle;
	for (const auto& tooth:teeth) {
		for (city_id i:tooth) {
			isPureHandle[i] = false;
		}
	}
	std::vector<std::pair<city_id, city_id>> oneNeighbors = findOneNeighbors(isPureHandle, lpData, solution);
	simplifyPureOnePaths(oneNeighbors, isHandle);
	//Knoten, die beim Vereinfachen aus dem Griff entfernt wurden, auch tatsächlich entfernen
	for (size_t i = 0; i < handle.size();) {
		if (isHandle[handle[i]]) {
			++i;
		} else {
			handle[i] = handle.back();
			handle.pop_back();
		}
	}
	//Prüfen, dass der Kamm noch gültig ist
	validate(cityCount);
}

/*
 * Findet alle 1-Kanten im reinen Griff und gibt sie in folgender Form zurück: Falls ein Knoten nicht im reinen Griff
 * liegt oder keine 1-Nachbarn im reinen Griff hat, sind beide Werte TSPInstance::invalid_city. Falls er genau einen
 * solchen Nachbarn hat, ist dies der erste Wert, der zweite ist invalid_city. Sonst sind die Werte die beiden Nachbarn.
 */
std::vector<std::pair<city_id, city_id>> CombHeuristic::Comb::findOneNeighbors(const std::vector<bool>& isPureHandle,
																			   const TspLpData& lpData,
																			   const std::vector<double>& solution) {
	const city_id invalid_city = TSPInstance::invalid_city;
	city_id cityCount = lpData.getTSP().getCityCount();
	std::vector<std::pair<city_id, city_id>> oneNeighbors(cityCount, {invalid_city, invalid_city});
	lemon::Tolerance<double> localTolerance;
	for (city_id i:handle) {
		if (isPureHandle[i]) {
			std::vector<city_id> validNeighbors;
			for (city_id neighbor = 0; neighbor < cityCount; ++neighbor) {
				if (isPureHandle[neighbor]) {
					variable_id var = lpData.getVariable(i, neighbor);
					if (var != LinearProgram::invalid_variable) {
						if (!localTolerance.different(solution[var], 1)) {
							validNeighbors.push_back(neighbor);
						}
					}
				}
			}
			if (validNeighbors.size() == 2) {
				oneNeighbors[i] = {validNeighbors[0], validNeighbors[1]};
			} else if (validNeighbors.size() == 1) {
				oneNeighbors[i] = {validNeighbors[0], TSPInstance::invalid_city};
			}
		}
	}
	return oneNeighbors;
}

/**
 * Vereinfacht 1-Pfade mit mindestens 4 Knoten im reinen Griff durch das Erstellen von zwei neuen Blättern, siehe
 * Grötschel, Holland 1991, Seite 25 bzw 165. Die neuen Zähne werden hinzugefügt, der Griff wird noch nicht verändert.
 * @param oneNeighbors der von findOneNeighbors erzeugt Vector von 1-Nachbarn im reinen Griff. Wird überschrieben.
 * @param isHandle gibt sowohl am Anfang als auch am Ende des Aufrufs an, welche Knoten im Griff liegen (sollten)
 */
void CombHeuristic::Comb::simplifyPureOnePaths(std::vector<std::pair<city_id, city_id>>& oneNeighbors,
											   std::vector<bool>& isHandle) {
	const city_id invalid_city = TSPInstance::invalid_city;
	for (city_id pathStart:handle) {
		std::pair<city_id, city_id>& startNeighbors = oneNeighbors[pathStart];
		//Der Knoten hat genau einen 1-Nachbarn im reinen Griff, kann also Anfang eines 1-Pfades sein
		if (startNeighbors.first != invalid_city && startNeighbors.second == invalid_city) {
			if (oneNeighbors[startNeighbors.first].second == invalid_city) {
				//Der 1-Nachbar hat auch nur einen 1-Nachbarn (zwangsläufig pathStart), der Pfad hat nur Länge 1
				//Also kann keiner der beiden Knoten in einem Pfad genutzt werden
				oneNeighbors[startNeighbors.first].first = invalid_city;
				startNeighbors.first = invalid_city;
				continue;
			}
			city_id pathEnd = pathStart;
			city_id nextInPath = startNeighbors.first;
			while (nextInPath != invalid_city) {
				if (pathEnd != pathStart) {
					//Es gibt einen nächsten Knoten, also ist pathEnd ein innerer Knoten und wird aus dem Griff entfernt
					//Außerdem kann pathEnd nicht in weiteren Pfaden vorkommen
					isHandle[pathEnd] = false;
					oneNeighbors[pathEnd] = {invalid_city, invalid_city};
				}
				city_id next;
				std::pair<city_id, city_id> neighborsCurr = oneNeighbors[nextInPath];
				if (neighborsCurr.first != pathEnd) {
					next = neighborsCurr.first;
				} else {
					next = neighborsCurr.second;
				}
				pathEnd = nextInPath;
				nextInPath = next;
				assert(nextInPath != pathStart);
			}
			assert(pathEnd != pathStart);
			//Der zweite und der vorletzte Knoten im Pfad
			city_id secondNode = startNeighbors.first;
			city_id secondToLast = oneNeighbors[pathEnd].first;
			if (secondNode == secondToLast) {
				//Falls sie gleich sind, hat der Pfad nur 3 Knoten und kann nicht vereinfacht werden
				isHandle[secondNode] = true;
			} else {
				//Falls sie ungleich sind, kann der Kamm vereinfacht werden
				assert(isHandle[pathStart]);
				assert(isHandle[pathEnd]);
				assert(!isHandle[secondNode]);
				assert(!isHandle[secondToLast]);
				teeth.push_back({pathStart, secondNode});
				teeth.push_back({pathEnd, secondToLast});
			}
			//Anfang und Ende können nicht weiter genutzt werden
			startNeighbors = oneNeighbors[pathEnd] = {invalid_city, invalid_city};
		}
	}
}

//Prüft, dass der Kamm "gültig" (d.h. tatsächlich ein Kamm) ist
void CombHeuristic::Comb::validate(city_id cityCount) const {
	//Falls NDEBUG definiert ist, ist assert "deaktiviert"
#ifndef NDEBUG
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
			inComb[i] = true;
		}
		assert(foundHandle);
		assert(foundNonHandle);
	}
#endif
}
