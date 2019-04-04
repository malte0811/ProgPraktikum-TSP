#include <connectivity_cut_gen.hpp>
#include <stack>
#include <tsp_lp_data.hpp>

ConnectivityCutGen::ConnectivityCutGen(const TSPInstance& inst, const TspLpData& lpData)
		: tsp(inst), workToOrig(workGraph), lpData(lpData), origToWork(inst.getCityCount()) {
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

CutGenerator::CutStatus ConnectivityCutGen::validate(LinearProgram& lp, const std::vector<double>& solution,
													 CutGenerator::CutStatus currentStatus) {
	//Grundzustand wiederherstellen
	for (Graph::EdgeIt it(workGraph); it != lemon::INVALID;) {
		Graph::Edge e = it;
		++it;
		workGraph.erase(e);
	}
	/*
	 * Alle Kanten einfügen, deren Variablen einen echt positiven Wert haben (Kanten mit Wert 0 ändern den minimalen
	 * Schnitt nicht)
	 */
	for (variable_id i = 0; i < solution.size(); ++i) {
		if (tolerance.positive(solution[i])) {
			TspLpData::Edge e = lpData.getEdge(i);
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
				std::vector<variable_id> indices;
				const std::vector<city_id>& currentComponent = components[i];
				for (size_t aId = 1; aId < currentComponent.size(); ++aId) {
					for (size_t bId = 0; bId < aId; ++bId) {
						variable_id var = lpData.getVariable(currentComponent[aId], currentComponent[bId]);
						if (var != LinearProgram::invalid_variable) {
							indices.push_back(var);
						}
					}
				}
				assert(!indices.empty());
				constrs.emplace_back(indices, std::vector<double>(indices.size(), 1), LinearProgram::less_eq,
									 currentComponent.size() - 1);
			}
		}
		lp.addConstraints(constrs);
		return recalc;
	} else {
		return valid;
	}
}
