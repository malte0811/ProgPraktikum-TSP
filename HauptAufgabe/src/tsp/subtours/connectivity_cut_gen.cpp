#include <connectivity_cut_gen.hpp>
#include <lemon/core.h>
#include <cstddef>
#include <algorithm>
#include <stack>
#include <tsp_lp_data.hpp>
#include <cut_generator.hpp>
#include <linear_program.hpp>
#include <tsp_instance.hpp>

using std::size_t;

ConnectivityCutGen::ConnectivityCutGen(const TSPInstance& inst, const TspLpData& lpData)
		: origToWork(inst.getCityCount()), workToOrig(workGraph), tsp(inst), lpData(lpData) {
	/*
	 * Grundzustand des Arbeitsgraphen herstellen: So viele Knoten wie in der TSP-Instanz, aber keine Kanten.
	 * Außerdem NodeMap's anlegen, um Städten in der TSP-Instanten in Knoten im Arbeitsgraphen umzuwandeln und umgekehrt
	 */
	for (city_id i = 0; i < tsp.getCityCount(); ++i) {
		Graph::Node newNode = workGraph.addNode();
		origToWork[i] = newNode;
		workToOrig[newNode] = i;
	}
}

CutGenerator::CutStatus ConnectivityCutGen::validate(LinearProgram& lp, const std::vector<double>& solution,
													 CutGenerator::CutStatus) {
	//Grundzustand wiederherstellen
	for (Graph::EdgeIt it(workGraph); it != lemon::INVALID;) {
		Graph::Edge e = it;
		++it;
		workGraph.erase(e);
	}
	/*
	 * Alle Kanten einfügen, deren Variablen einen echt positiven Wert haben
	 */
	for (variable_id i = 0; i < lpData.getVariableCount(); ++i) {
		if (tolerance.positive(solution[i])) {
			const TspLpData::Edge& e = lpData.getEdge(i);
			workGraph.addEdge(origToWork[e.first], origToWork[e.second]);
		}
	}
	Graph::NodeMap<bool> visited(workGraph);
	std::vector<std::vector<city_id>> components;
	size_t maxIndex = 0;
	size_t maxSize = 0;
	//Graphendurchmusterung, um die Zusammenhangskomponenten zu bestimmen
	for (Graph::NodeIt startIt(workGraph); startIt != lemon::INVALID; ++startIt) {
		if (!visited[startIt]) {
			std::stack<Graph::Node> open;
			std::vector<city_id> currentComponent{workToOrig[startIt]};
			open.push(startIt);
			visited[startIt] = true;
			while (!open.empty()) {
				Graph::Node current = open.top();
				open.pop();
				for (Graph::OutArcIt nextIt(workGraph, current); nextIt != lemon::INVALID; ++nextIt) {
					Graph::Node neighbor = workGraph.target(nextIt);
					if (!visited[neighbor]) {
						visited[neighbor] = true;
						open.push(neighbor);
						currentComponent.push_back(workToOrig[neighbor]);
					}
				}
			}
			//Komponente mit den meisten Knoten speichern
			if (currentComponent.size() > maxSize) {
				maxSize = currentComponent.size();
				maxIndex = components.size();
			}
			components.push_back(currentComponent);
		}
	}
	if (components.size() > 1) {
		std::vector<LinearProgram::Constraint> constrs;
		constrs.reserve(components.size() - 1);
		/*
		 * Constraints für alle Komponenten außer der größten hinzufügen: Die letzte Constraint wird von den anderen
		 * impliziert, muss also nicht hinzugefügt werden. Außerdem ist der LP-Solver schneller, je dünner die Constraints
		 * besetzt sind (je größer die Zusammenhangskomponente, desto dichter die Constraint)
		 */
		for (size_t i = 0; i < components.size(); ++i) {
			if (i != maxIndex) {
				const std::vector<city_id>& currentComponent = components[i];
				std::vector<variable_id> usedVars;
				std::vector<double> coeffs(lpData.getVariableCount());
				city_id rhs = lpData.sparserInducedSum(currentComponent, coeffs, usedVars);
				constrs.push_back(LinearProgram::Constraint::fromDense(usedVars, coeffs, LinearProgram::less_eq, rhs));
			}
		}
		lp.addConstraints(constrs.begin(), constrs.end());
		return recalc;
	} else {
		return valid;
	}
}
