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
#include <comb_cut_gen.hpp>
#include <connectivity_cut_gen.hpp>
#include <tsp_lp_data.hpp>

namespace tspsolvers {

	/**
	 * Berechnet eine "relative kurze" Tour in der gegebenen TSP-Instanz, indem immer die kürzeste mögliche Kante
	 * hinzugefügt wird.
	 */
	TSPSolution solveGreedy(const TSPInstance& inst) {
		TspLpData data(inst);//TODO zu TSPInstance Methoden für IDs (für allg. Verwendung) hinzufügen
		std::vector<variable_id> sortedEdges;
		for (variable_id i = 0; i<inst.getEdgeCount(); ++i) {
			sortedEdges.push_back(i);
		}
		//Sortieren nach Kosten der entsprechenden Kanten
		std::sort(sortedEdges.begin(), sortedEdges.end(),
				  [&data](variable_id edgeA, variable_id edgeB) {
					  return data.getCost(edgeA) < data.getCost(edgeB);
				  }
		);
		/*
		 * Falls an Knoten i weniger als 2 Kanten anliegen, gibt otherEnd[i] den anderen Knoten in der selben
		 * Zusammenhangskomponente wie i an, der auch Grad <2 hat (Falls i Grad 0 hat, ist dies i selbst), d.h. das
		 * andere Ende des Toursegments. Falls an i 2 Kanten anliegen, ist der Wert beliebig.
		 */
		std::vector<city_id> otherEnd(inst.getCityCount());
		for (city_id i = 0; i<inst.getCityCount(); ++i) {
			otherEnd[i] = i;
		}
		std::vector<bool> used(static_cast<size_t>(inst.getEdgeCount()));
		std::vector<size_t> degree(inst.getCityCount());
		unsigned addedEdges = 0;
		for (variable_id eId:sortedEdges) {
			//Falls an einem der beiden Enden schon 2 Kanten anliegen, kann die Kante nicht hinzugefügt werden
			TspLpData::Edge e = data.getEdge(eId);
			size_t& edgesAtU = degree[e.first];
			if (edgesAtU>=2) {
				continue;
			}
			size_t& edgesAtV = degree[e.second];
			if (edgesAtV>=2) {
				continue;
			}
			/*
			 * Wenn die Enden beide Grad <2 haben (also Enden von Toursegmenten sind), können sie genau dann verbunden
			 * werden, wenn sie nicht zum selben Segment gehören.
			 */
			if (otherEnd[e.first] != e.second) {
				++addedEdges;
				++edgesAtU;
				++edgesAtV;
				used[eId] = true;
				if (addedEdges==inst.getCityCount()-1) {
					break;//Tour ist fast vollständig, die letzte Kante ist aber eindeutig bestimmt
				}
				//Enden des neuen Segments setzen
				//newEndU/V zwischenspeichern, falls die Segmente aus einzelnen Knoten bestehen
				city_id newEndU = otherEnd[e.first];
				city_id newEndV = otherEnd[e.second];
				otherEnd[newEndV] = newEndU;
				otherEnd[newEndU] = newEndV;
			}
		}
		closeHamiltonPath(inst, used, degree, data);
		return TSPSolution(inst, used, data);
	}

	/**
	 * Fügt die fehlende Kante in einen Hamilton-Pfad in der TSP-Instanz inst ein
	 */
	void closeHamiltonPath(const TSPInstance& instance, std::vector<bool>& used, const std::vector<size_t>& degree,
			//TODO weg
						   const TspLpData& data) {
		city_id firstEnd = TSPInstance::invalid_city;
		for (city_id currCity = 0; currCity<instance.getCityCount(); ++currCity) {
			if (degree[currCity]==1) {
				if (firstEnd!=TSPInstance::invalid_city) {
					used[data.getVariable(firstEnd, currCity)] = true;
					break;
				} else {
					firstEnd = currCity;
				}
			}
		}
	}

	/**
	 * Löst eine TSP-Instanz exakt durch das Lösen eines linearen ganzzahligen Programms
	 * @param inst Die zu lösende TSP-Instanz
	 * @param initial Eine erste obere Schranke, oder nullptr falls ohne eine solche gearbeitet werden soll
	 * @param lpEnv Die zu verwendende LP-Umgebung
	 * @param maxOpenSize Die maximale Größe der offenen Menge in Bytes. Auf 0 setzen, um DFS zu verwenden
	 * @return eine optimal Lösung der gegebenen TSP-Instanz
	 */
	TSPSolution solveLP(const TSPInstance& inst, const TSPSolution* initial, CPXENVptr& lpEnv, size_t maxOpenSize) {
		TspLpData data(inst);
		data.setupLowerBounds();
		ConnectivityCutGen connected(inst, data);
		SubtourCutGen subtours(inst, data);
		TwoMatchingCutGen matchings(inst, data, true);
		CombCutGen combs(inst, data);
		LinearProgram lp(lpEnv, inst.getName(), LinearProgram::minimize);
		data.setupBasicLP(lp);
		BranchAndCut bac(lp, {&connected, &subtours, &combs, &matchings}, &data, maxOpenSize);
		if (initial!=nullptr) {
			std::vector<long> asVars(static_cast<size_t>(inst.getEdgeCount()));
			for (variable_id i = 0; i<inst.getEdgeCount(); ++i) {
				asVars[i] = initial->getVariable(i);
			}
			bac.setUpperBound(asVars, initial->getCost());
		}
		clock_t start = std::clock();
		bac.solve();
		clock_t end = std::clock();
		double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
		std::cout << "Branch and cut took " << elapsed_secs << " seconds" << std::endl;
		return data.getUpperBound();
	}
}