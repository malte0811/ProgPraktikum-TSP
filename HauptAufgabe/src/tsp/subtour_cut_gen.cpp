#include <subtour_cut_gen.hpp>
#include <lemon_fixes/nagamochi_ibaraki.h>
#include <cmath>
#include <stack>

SubtourCutGen::SubtourCutGen(const TSPInstance& inst)
		: tsp(inst), origToWork(static_cast<size_t>(tsp.getCityCount())), workToOrig(workGraph), capacity(workGraph),
		  minCut(workGraph, capacity), tolerance(1e-5) {
	/*
	 * Grundzustand des Arbeitsgraphen abspeichern: So viele Knoten wie in der TSP-Instanz, aber keine Kanten.
	 * Außerdem NodeMap's anlegen, um Knoten im TSP-Graphen in Knoten im Arbeitsgraphen umzuwandeln und umgekehrt
	 */
	for (city_id i = 0; i<tsp.getCityCount(); ++i) {
		Graph::Node newNode = workGraph.addNode();
		origToWork[i] = newNode;
		workToOrig[newNode] = i;
	}
	baseState.save(workGraph);
}

CutGenerator::CutStatus SubtourCutGen::validate(LinearProgram& lp, const std::vector<double>& solution,
												CutStatus currentStatus) {
	if (currentStatus != valid) {
		return valid;
	}
	//Grundzustand wiederherstellen und speichern
	baseState.restore();
	baseState.save(workGraph);
	/*
	 * Alle Kanten einfügen, deren Variablen einen echt positiven Wert haben (Kanten mit Wert 0 ändern den minimalen
	 * Schnitt nicht)
	 */
	for (city_id lowerEnd = 0; lowerEnd < tsp.getCityCount() - 1; ++lowerEnd) {
		for (city_id higherEnd = lowerEnd + 1; higherEnd < tsp.getCityCount(); ++higherEnd) {
			variable_id edgeVar = tsp.getVariable(higherEnd, lowerEnd);
			if (tolerance.positive(solution[edgeVar])) {
				Graph::Edge inWork = workGraph.addEdge(origToWork[lowerEnd], origToWork[higherEnd]);
				capacity[inWork] = solution[edgeVar];
			}
		}
	}
	minCut.run();
	double cutCapacity = minCut.minCutValue();
	if (!tolerance.less(cutCapacity, 2)) {
		return valid;
	}
	//Positive Kapazität <2->Verletzte Subtour-Constraint
	Graph::NodeMap<bool> inCut(workGraph);
	minCut.minCutMap(inCut);
	city_id cutSize = 0;
	for (Graph::NodeIt it(workGraph); it != lemon::INVALID; ++it) {
		if (inCut[it]) {
			++cutSize;
		}
	}
	/*
	 * Falls der gefundene Cut mehr als die Hälfte aller Knoten enthält, wird das Komplement hinzugefügt. Wie oben
	 * wird so die Anzahl an Einträgen ungleich 0 minimiert.
	 */
	const bool cutVal = cutSize < tsp.getCityCount() / 2;
	if (!cutVal) {
		cutSize = tsp.getCityCount() - cutSize;
	}
	//Kann auftreten, wenn alle Constraints erfüllt sind, aber Kanten mit sehr kleinen positiven Werten existieren
	if (cutSize <= 2) {
		return valid;
	}
	std::vector<int> induced;
	for (city_id lowerEnd = 0; lowerEnd < tsp.getCityCount() - 1; ++lowerEnd) {
		for (city_id higherEnd = lowerEnd + 1; higherEnd < tsp.getCityCount(); ++higherEnd) {
			bool vInCut = inCut[origToWork[lowerEnd]];
			bool uInCut = inCut[origToWork[higherEnd]];
			if (vInCut == cutVal && uInCut == cutVal) {
				induced.push_back(tsp.getVariable(higherEnd, lowerEnd));
			}
		}
	}
	assert(!induced.empty());
	lp.addConstraint(
			LinearProgram::Constraint(induced, std::vector<double>(induced.size(), 1), LinearProgram::less_eq,
									  cutSize - 1));
	return recalc;
}