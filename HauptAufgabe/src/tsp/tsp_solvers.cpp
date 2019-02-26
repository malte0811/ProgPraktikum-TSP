#include <iostream>
#include <ctime>
#include <istream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <stack>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <tsp_solution.hpp>
#include <tsp_solvers.hpp>
#include <subtour_cut_gen.hpp>
#include <two_matching_cut_gen.hpp>
#include <branch_and_cut.hpp>

namespace tspsolvers {

	/**
		 * Berechnet eine "relative kurze" Tour in der gegebenen TSP-Instanz, indem immer die kürzeste mögliche Kante
		 * hinzugefügt wird.
		 */
	TSPSolution solveGreedy(const TSPInstance& inst) {
		const Graph::EdgeMap <cost_t>& distances = inst.getGraphDistances();
		std::vector<Graph::Edge> sortedEdges;
		for (Graph::EdgeIt it(inst.getGraph()); it!=lemon::INVALID; ++it) {
			sortedEdges.push_back(it);
		}
		//Sortieren nach Kosten der entsprechenden Kanten
		std::sort(sortedEdges.begin(), sortedEdges.end(),
				  [&distances](const Graph::Edge& edgeA, const Graph::Edge& edgeB) {
					  return distances[edgeA]<distances[edgeB];
				  }
		);
		/*
		 * Falls an Knoten i weniger als 2 Kanten anliegen, gibt otherEnd[i] den anderen Knoten in der selben
		 * Zusammenhangskomponente wie i an, der auch Grad <2 hat (Falls i Grad 0 hat, ist dies i selbst), d.h. das
		 * andere Ende des Toursegments. Falls an i 2 Kanten anliegen, ist der Wert beliebig.
		 */
		Graph::NodeMap <Graph::Node> otherEnd(inst.getGraph());
		for (Graph::NodeIt it(inst.getGraph()); it!=lemon::INVALID; ++it) {
			otherEnd[it] = it;
		}
		std::vector<bool> used(static_cast<size_t>(inst.getEdgeCount()));
		Graph::NodeMap <size_t> degree(inst.getGraph());
		unsigned addedEdges = 0;
		for (const Graph::Edge& e:sortedEdges) {
			//Falls an einem der beiden Enden schon 2 Kanten anliegen, kann die Kante nicht hinzugefügt werden
			Graph::Node u = inst.getGraph().u(e);
			size_t& edgesAtU = degree[u];
			if (edgesAtU>=2) {
				continue;
			}
			Graph::Node v = inst.getGraph().v(e);
			size_t& edgesAtV = degree[v];
			if (edgesAtV>=2) {
				continue;
			}
			/*
			 * Wenn die Enden beide Grad <2 haben (also Enden von Toursegmenten sind), können sie genau dann verbunden
			 * werden, wenn sie nicht zum selben Segment gehören.
			 */
			if (otherEnd[u]!=v) {
				++addedEdges;
				++edgesAtU;
				++edgesAtV;
				used[inst.getVariable(e)] = true;
				if (addedEdges==inst.getSize()-1) {
					break;//Tour ist fast vollständig, die letzte Kante ist aber eindeutig bestimmt
				}
				//Enden des neuen Segments setzen
				//newEndU/V zwischenspeichern, falls die Segmente aus einzelnen Knoten bestehen
				Graph::Node newEndU = otherEnd[u];
				Graph::Node newEndV = otherEnd[v];
				otherEnd[newEndV] = newEndU;
				otherEnd[newEndU] = newEndV;
			}
		}
		closeHamiltonPath(inst, used, degree);
		return TSPSolution(inst, used);
	}

	/**
	 * Fügt die fehlende Kante in einen Hamilton-Pfad in der TSP-Instanz inst ein
	 */
	void closeHamiltonPath(const TSPInstance& instance, std::vector<bool>& used,
						   const Graph::NodeMap <size_t>& degree) {
		Graph::Node firstEnd = lemon::INVALID;
		for (Graph::NodeIt it(instance.getGraph()); it!=lemon::INVALID; ++it) {
			if (degree[it]==1) {
				if (firstEnd!=lemon::INVALID) {
					city_id cityA = instance.getCity(firstEnd);
					city_id cityB = instance.getCity(it);
					used[instance.getVariable(cityA, cityB)] = true;
					break;
				} else {
					firstEnd = it;
				}
			}
		}
	}

	TSPSolution solveLP(const TSPInstance& inst, const TSPSolution* initial, CPXENVptr& lpEnv, size_t maxOpenSize) {
		LinearProgram lp(lpEnv, inst.getName(), LinearProgram::minimize);
		inst.setupBasicLP(lp);
		SubtourCutGen subtours(inst);
		TwoMatchingCutGen matchings(inst, true);
		BranchAndCut bac(lp, {&subtours, &matchings}, maxOpenSize);
		if (initial!=nullptr) {
			std::vector<long> asVars(static_cast<size_t>(inst.getEdgeCount()));
			for (variable_id i = 0; i<inst.getEdgeCount(); ++i) {
				asVars[i] = initial->getVariable(i);
			}
			bac.setUpperBound(asVars, initial->getCost());
		}
		std::vector<long> tour = bac.solve();
		std::vector<bool> asBools(tour.size());
		for (size_t i = 0; i<tour.size(); ++i) {
			asBools[i] = tour[i];
		}
		return TSPSolution(inst, asBools);
	}
}