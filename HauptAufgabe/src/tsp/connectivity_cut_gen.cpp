#include <connectivity_cut_gen.hpp>
#include <stack>

ConnectivityCutGen::ConnectivityCutGen(const TSPInstance& inst) : tsp(inst), workToOrig(workGraph),
																  origToWork(inst.getCityCount()) {
	/*
	 * Grundzustand des Arbeitsgraphen abspeichern: So viele Knoten wie in der TSP-Instanz, aber keine Kanten.
	 * Außerdem NodeMap's anlegen, um Knoten im TSP-Graphen in Knoten im Arbeitsgraphen umzuwandeln und umgekehrt
	 */
	for (city_id i = 0; i < tsp.getCityCount(); ++i) {
		Graph::Node newNode = workGraph.addNode();
		origToWork[i] = newNode;
		workToOrig[newNode] = i;
	}
	baseState.save(workGraph);
}

CutGenerator::CutStatus ConnectivityCutGen::validate(LinearProgram& lp, const std::vector<double>& solution,
													 CutGenerator::CutStatus currentStatus) {
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
				workGraph.addEdge(origToWork[lowerEnd], origToWork[higherEnd]);
			}
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
						indices.push_back(tsp.getVariable(currentComponent[aId], currentComponent[bId]));
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
