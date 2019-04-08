#include <subtour_cut_gen.hpp>
#include <lemon_fixes/nagamochi_ibaraki.h>
#include <cmath>
#include <stack>
#include <tsp_lp_data.hpp>
#include <tsp_utils.hpp>

SubtourCutGen::SubtourCutGen(const TSPInstance& inst, const TspLpData& lpData)
		: tsp(inst), origToWork(static_cast<size_t>(tsp.getCityCount())), workToOrig(workGraph), capacity(workGraph),
		  minCut(workGraph, capacity), tolerance(1e-5), lpData(lpData) {
	/*
	 * Grundzustand des Arbeitsgraphen abspeichern: So viele Knoten wie in der TSP-Instanz, aber keine Kanten.
	 * Außerdem NodeMap's anlegen, um Knoten im TSP-Graphen in Knoten im Arbeitsgraphen umzuwandeln und umgekehrt
	 */
	for (city_id i = 0; i<tsp.getCityCount(); ++i) {
		Graph::Node newNode = workGraph.addNode();
		origToWork[i] = newNode;
		workToOrig[newNode] = i;
	}
}

CutGenerator::CutStatus SubtourCutGen::validate(LinearProgram& lp, const std::vector<double>& solution,
												CutStatus currentStatus) {
	if (currentStatus != valid) {
		return valid;
	}
	tsp_util::addSupportGraphEdges(tsp, lpData, tolerance, solution, workGraph, origToWork, workToOrig, capacity);
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
	std::vector<variable_id> induced;
	for (city_id lowerEnd = 0; lowerEnd < tsp.getCityCount() - 1; ++lowerEnd) {
		for (city_id higherEnd = lowerEnd + 1; higherEnd < tsp.getCityCount(); ++higherEnd) {
			bool vInCut = inCut[origToWork[lowerEnd]];
			bool uInCut = inCut[origToWork[higherEnd]];
			if (vInCut == cutVal && uInCut == cutVal) {
				variable_id var = lpData.getVariable(higherEnd, lowerEnd);
				if (var >= 0) {
					induced.push_back(var);
				}
			}
		}
	}
	assert(!induced.empty());
	lp.addConstraint(
			LinearProgram::Constraint(induced, std::vector<double>(induced.size(), 1), LinearProgram::less_eq,
									  cutSize - 1));
	return recalc;
}