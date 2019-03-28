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
	//Grundzustand wiederherstellen und speichern
	baseState.restore();
	baseState.save(workGraph);
	/*
	 * Alle Kanten einfügen, deren Variablen einen echt positiven Wert haben (Kanten mit Wert 0 ändern den minimalen
	 * Schnitt nicht)
	 */
	for (city_id lowerEnd = 0; lowerEnd<tsp.getCityCount()-1; ++lowerEnd) {
		for (city_id higherEnd = lowerEnd+1; higherEnd<tsp.getCityCount(); ++higherEnd) {
			variable_id edgeVar = tsp.getVariable(higherEnd, lowerEnd);
			if (tolerance.positive(solution[edgeVar])) {
				Graph::Edge inWork = workGraph.addEdge(origToWork[lowerEnd], origToWork[higherEnd]);
				capacity[inWork] = solution[edgeVar];
			}
		}
	}
	minCut.run();
	double capacity = minCut.minCutValue();
	if (capacity==0) {
		/*
		 * Da es keine Kanten mit Kapazität 0 gibt, kann der Graph nicht zusammenhängend sein, es können also leicht
		 * mehrere verletzte Constraints bestimmt werden.
		 */
		addConnectivityConstraints(lp);
		return CutGenerator::recalc;
	} else if (tolerance.less(capacity, 2)) {
		//Positive Kapazität <2->Verletzte Subtour-Constraint
		addCutConstraint(lp);
		return CutGenerator::recalc;
	} else {
		return CutGenerator::valid;
	}
}

/**
 * Fügt die Constraints hinzu, die den Zusammenhangskomponenten de Arbeitsgraphen entsprechen
 */
void SubtourCutGen::addConnectivityConstraints(LinearProgram& lp) {
	Graph::NodeMap<bool> visited(workGraph);
	std::vector<std::vector<city_id>> components;
	size_t maxIndex = 0;
	size_t maxSize = 0;
	//Graphendurchmusterung, um die Zusammenhangskomponenten zu bestimmen
	for (Graph::NodeIt startIt(workGraph); startIt!=lemon::INVALID; ++startIt) {
		if (!visited[startIt]) {
			std::stack<Graph::Node> open;
			std::vector<city_id> currentComponent{workToOrig[startIt]};
			open.push(startIt);
			visited[startIt] = true;
			while (!open.empty()) {
				Graph::Node current = open.top();
				open.pop();
				for (Graph::OutArcIt nextIt(workGraph, current); nextIt!=lemon::INVALID; ++nextIt) {
					Graph::Node neighbor = workGraph.target(nextIt);
					if (!visited[neighbor]) {
						visited[neighbor] = true;
						open.push(neighbor);
						currentComponent.push_back(workToOrig[neighbor]);
					}
				}
			}
			//Komponente mit den meisten Knoten speichern
			if (currentComponent.size()>maxSize) {
				maxSize = currentComponent.size();
				maxIndex = components.size();
			}
			components.push_back(currentComponent);
		}
	}
	std::vector<LinearProgram::Constraint> constrs;
	constrs.reserve(components.size()-1);
	/*
	 * Constraints für alle Komponenten außer der größten hinzufügen: Die letzte Constraint wird von den anderen
	 * impliziert, muss also nicht hinzugefügt werden. Außerdem ist der LP-Solver schneller, je dünner die Constraints
	 * besetzt sind (je größer die Zusammenhangskomponente, desto dichter die Constraint)
	 */
	for (size_t i = 0; i<components.size(); ++i) {
		if (i!=maxIndex) {
			std::vector<variable_id> indices;
			const std::vector<city_id>& currentComponent = components[i];
			for (size_t aId = 1; aId<currentComponent.size(); ++aId) {
				for (size_t bId = 0; bId<aId; ++bId) {
					indices.push_back(tsp.getVariable(currentComponent[aId], currentComponent[bId]));
				}
			}
			constrs.emplace_back(indices, std::vector<double>(indices.size(), 1), LinearProgram::less_eq,
								 currentComponent.size()-1);
		}
	}
	lp.addConstraints(constrs);
}

/**
 * Fügt die Constraint hinzu, die dem (bereits berechneten) Min-Cut entspricht
 */
void SubtourCutGen::addCutConstraint(LinearProgram& lp) {
	Graph::NodeMap<bool> inCut(workGraph);
	minCut.minCutMap(inCut);
	city_id cutSize = 0;
	for (Graph::NodeIt it(workGraph); it!=lemon::INVALID; ++it) {
		if (inCut[it]) {
			++cutSize;
		}
	}
	/*
	 * Falls der gefundene Cut mehr als die Hälfte aller Knoten enthält, wird das Komplement hinzugefügt. Wie oben
	 * wird so die Anzahl an Einträgen ungleich 0 minimiert.
	 */
	const bool cutVal = cutSize<tsp.getCityCount()/2;
	if (!cutVal) {
		cutSize = tsp.getCityCount()-cutSize;
	}
	std::vector<int> induced;
	for (city_id lowerEnd = 0; lowerEnd<tsp.getCityCount()-1; ++lowerEnd) {
		for (city_id higherEnd = lowerEnd+1; higherEnd<tsp.getCityCount(); ++higherEnd) {
			bool vInCut = inCut[origToWork[lowerEnd]];
			bool uInCut = inCut[origToWork[higherEnd]];
			if (vInCut==cutVal && uInCut==cutVal) {
				induced.push_back(tsp.getVariable(higherEnd, lowerEnd));
			}
		}
	}
	lp.addConstraint(LinearProgram::Constraint(induced, std::vector<double>(induced.size(), 1), LinearProgram::less_eq,
											   cutSize-1));
}