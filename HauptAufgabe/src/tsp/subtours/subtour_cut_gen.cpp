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
	if (currentStatus == recalc) {
		//Der Zusammenhangs-CutGen hat schon verletzte Subtour-Ungleichungen gefunden
		return valid;
	}
	tsp_util::addSupportGraphEdges(tsp, lpData, tolerance, solution, workGraph, origToWork, workToOrig, capacity);
	minCut.run();
	double cutCapacity = minCut.minCutValue();
	const double cutScale = 100;
	if (!tolerance.less(cutScale * cutCapacity, cutScale * 2)) {
		return valid;
	}
	//Positive Kapazität<2 ->Verletzte Subtour-Constraint
	Graph::NodeMap<bool> inCut(workGraph);
	minCut.minCutMap(inCut);
	std::vector<city_id> inducing;
	for (Graph::NodeIt it(workGraph); it != lemon::INVALID; ++it) {
		if (inCut[it]) {
			inducing.push_back(workToOrig[it]);
		}
	}
	//Kann auftreten, wenn alle Constraints erfüllt sind, aber Kanten mit sehr kleinen positiven Werten existieren
	if (inducing.size() <= 2 || inducing.size() >= tsp.getCityCount() - 2) {
		return valid;
	}
	std::vector<variable_id> usedVars;
	std::vector<double> coeffs(lpData.getVariableCount());
	city_id rhs = lpData.sparserInducedSum(inducing, coeffs, usedVars);
	LinearProgram::Constraint constr = LinearProgram::Constraint::fromDense(usedVars, coeffs, LinearProgram::less_eq,
																			rhs);
	//Möglich, da Kanten mit kleinen positiven Werte nicht zum Graphen hinzugefügt werden
	if (constr.isViolated(solution, tolerance)) {
		lp.addConstraint(constr);
		return recalc;
	} else {
		return valid;
	}
}