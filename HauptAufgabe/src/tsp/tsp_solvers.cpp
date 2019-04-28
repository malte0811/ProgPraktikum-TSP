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
#include <connectivity_cut_gen.hpp>
#include <tsp_lp_data.hpp>
#include <comb_cut_gen.hpp>

namespace tspsolvers {

	namespace cutgens {
		extern const char *const connected = "connected";
		extern const char *const subtour = "subtour";
		extern const char *const twoMatching = "2matching";
		extern const char *const generalCombs = "combs";
		extern const char *const defaultGens = "connected,subtour,2matching,combs";
	}

	/**
	 * Berechnet eine "relative kurze" Tour in der gegebenen TSP-Instanz, indem immer die kürzeste mögliche Kante
	 * hinzugefügt wird.
	 */
	TSPSolution solveGreedy(const TSPInstance& inst) {
		using lemon::FullGraph;
		FullGraph g(inst.getCityCount());
		std::vector<FullGraph::Edge> sortedEdges;
		for (FullGraph::EdgeIt it(g); it != lemon::INVALID; ++it) {
			sortedEdges.push_back(it);
		}
		//Sortieren nach Kosten der entsprechenden Kanten
		std::sort(sortedEdges.begin(), sortedEdges.end(),
				  [&inst, &g](FullGraph::Edge e1, FullGraph::Edge e2) {
					  city_id endA1 = FullGraph::id(g.u(e1));
					  city_id endA2 = FullGraph::id(g.u(e2));
					  city_id endB1 = FullGraph::id(g.v(e1));
					  city_id endB2 = FullGraph::id(g.v(e2));
					  return inst.getDistance(endA1, endB1) < inst.getDistance(endA2, endB2);
				  }
		);
		/*
		 * Falls an Knoten i weniger als 2 Kanten anliegen, gibt otherEnd[i] den anderen Knoten in der selben
		 * Zusammenhangskomponente wie i an, der auch Grad <2 hat (Falls i Grad 0 hat, ist dies i selbst), d.h. das
		 * andere Ende des Toursegments. Falls an i 2 Kanten anliegen, ist der Wert beliebig.
		 */
		FullGraph::NodeMap<FullGraph::Node> otherEnd(g);
		for (FullGraph::NodeIt it(g); it != lemon::INVALID; ++it) {
			otherEnd[it] = it;
		}
		FullGraph::EdgeMap<bool> used(g);
		FullGraph::NodeMap<size_t> degree(g);
		int addedEdges = 0;
		for (FullGraph::Edge e:sortedEdges) {
			//Falls an einem der beiden Enden schon 2 Kanten anliegen, kann die Kante nicht hinzugefügt werden
			size_t& edgesAtU = degree[g.u(e)];
			if (edgesAtU >= 2) {
				continue;
			}
			size_t& edgesAtV = degree[g.v(e)];
			if (edgesAtV >= 2) {
				continue;
			}
			/*
			 * Wenn die Enden beide Grad <2 haben (also Enden von Toursegmenten sind), können sie genau dann verbunden
			 * werden, wenn sie nicht zum selben Segment gehören.
			 */
			if (otherEnd[g.u(e)] != g.v(e)) {
				++addedEdges;
				++edgesAtU;
				++edgesAtV;
				used[e] = true;
				if (addedEdges == inst.getCityCount() - 1) {
					break;//Tour ist fast vollständig, die letzte Kante ist aber eindeutig bestimmt
				}
				//Enden des neuen Segments setzen
				//newEndU/V zwischenspeichern, falls die Segmente aus einzelnen Knoten bestehen
				FullGraph::Node newEndU = otherEnd[g.u(e)];
				FullGraph::Node newEndV = otherEnd[g.v(e)];
				otherEnd[newEndV] = newEndU;
				otherEnd[newEndU] = newEndV;
			}
		}
		closeHamiltonPath(g, used, degree);
		return TSPSolution(inst, g, used);
	}

	/**
	 * Fügt die fehlende Kante in einen Hamilton-Pfad in der TSP-Instanz inst ein
	 */
	void closeHamiltonPath(const lemon::FullGraph& g, lemon::FullGraph::EdgeMap<bool>& used,
						   const lemon::FullGraph::NodeMap<size_t>& degree) {
		using lemon::FullGraph;
		FullGraph::Node firstEnd = lemon::INVALID;
		for (FullGraph::NodeIt it(g); it != lemon::INVALID; ++it) {
			if (degree[it] == 1) {
				if (firstEnd != lemon::INVALID) {
					used[g.edge(firstEnd, it)] = true;
					break;
				} else {
					firstEnd = it;
				}
			}
		}
	}

	/**
	 * Löst eine TSP-Instanz exakt durch das Lösen eines linearen ganzzahligen Programms
	 * @param inst Die zu lösende TSP-Instanz
	 * @param initial Eine erste obere Schranke, oder nullptr falls ohne eine solche gearbeitet werden soll
	 * @param lpEnv Die zu verwendende LP-Umgebung
	 * @param dfs Gibt an, ob DFS verwendet werden soll
	 * @param cutGenerators die Namen der zu nutzenden CutGens (siehe tspsolvers::cutgens)
	 * @return eine optimal Lösung der gegebenen TSP-Instanz
	 */
	TSPSolution solveLP(const TSPInstance& inst, const TSPSolution *initial, CPXENVptr& lpEnv, bool dfs,
						const std::vector<std::string>& cutGenerators) {
		TspLpData data(inst, initial);
		std::cout << "Starting at variable count " << data.getVariableCount() << std::endl;
		ConnectivityCutGen connected(inst, data);
		SubtourCutGen subtours(inst, data);
		TwoMatchingCutGen matchings(inst, data, true);
		CombCutGen generalCombs(inst, data);
		LinearProgram lp(lpEnv, inst.getName(), LinearProgram::minimize);
		data.setupBasicLP(lp);
		std::vector<CutGenerator *> gens;
		gens.reserve(cutGenerators.size());
		for (const std::string& genName:cutGenerators) {
			if (tspsolvers::cutgens::connected == genName) {
				gens.push_back(&connected);
			} else if (tspsolvers::cutgens::subtour == genName) {
				gens.push_back(&subtours);
			} else if (tspsolvers::cutgens::twoMatching == genName) {
				gens.push_back(&matchings);
			} else if (tspsolvers::cutgens::generalCombs == genName) {
				gens.push_back(&generalCombs);
			} else {
				throw std::runtime_error("Unknown cut generator: " + genName);
			}
		}
		BranchAndCut bac(lp, gens, &data, dfs);
		if (initial != nullptr) {
			bac.setUpperBound({}, initial->getCost());
		}
		clock_t start = std::clock();
		bac.solve();
		clock_t end = std::clock();
		double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
		std::cout << "Branch and cut took " << elapsed_secs << " seconds" << std::endl;
		return data.getUpperBound();
	}
}