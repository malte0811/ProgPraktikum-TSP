#include <subtour_cut_gen.hpp>
#include <lemon/core.h>
#include <lemon_fixes/nagamochi_ibaraki.h>
#include <cstddef>
#include <tsp_lp_data.hpp>
#include <tsp_utils.hpp>
#include <cut_generator.hpp>
#include <linear_program.hpp>
#include <tsp_instance.hpp>

using std::size_t;

SubtourCutGen::SubtourCutGen(const TSPInstance& inst, const TspLpData& lpData)
		: tsp(inst), lpData(lpData), origToWork(static_cast<size_t>(tsp.getCityCount())), workToOrig(workGraph),
		  capacity(workGraph),
		  minCut(workGraph, capacity), tolerance(1e-5) {
	/*
	 * Grundzustand des Arbeitsgraphen abspeichern: So viele Knoten wie in der TSP-Instanz, aber keine Kanten.
	 * Außerdem NodeMap's anlegen, um Knoten im TSP-Graphen in Knoten im Arbeitsgraphen umzuwandeln und umgekehrt
	 */
	for (city_id i = 0; i < tsp.getCityCount(); ++i) {
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
	tsp_util::addSupportGraphEdges(lpData, tolerance, solution, workGraph, origToWork, workToOrig, capacity);
	minCut.run();
	double cutCapacity = minCut.minCutValue();
	if (!tolerance.less(cutCapacity, 2)) {
		return valid;
	}

	//Kapazität<2 ->Verletzte Subtour-Constraint
	Graph::NodeMap<bool> inCut(workGraph);
	minCut.minCutMap(inCut);
	std::vector<city_id> inducing;
	for (Graph::NodeIt it(workGraph); it != lemon::INVALID; ++it) {
		if (inCut[it]) {
			inducing.push_back(workToOrig[it]);
		}
	}

	//Kann auftreten, wenn alle Constraints erfüllt sind, aber Kanten mit sehr kleinen positiven Werten existieren
	if (inducing.size() <= 2 || static_cast<city_id>(inducing.size()) >= tsp.getCityCount() - 2) {
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