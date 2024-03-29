#include <comb_cut_gen.hpp>
#include <cassert>
#include <lemon/core.h>
#include <cstddef>
#include <array>
#include <comb_heuristic.hpp>
#include <iostream>
#include <set>
#include <contraction_rule.hpp>
#include <tsp_lp_data.hpp>
#include <tsp_utils.hpp>
#include <cut_generator.hpp>
#include <linear_program.hpp>
#include <tsp_instance.hpp>

using std::size_t;

CombCutGen::CombCutGen(const TSPInstance& tsp, const TspLpData& lpData) :
//ContractionRule's für die im Paper beschriebenen Transformationen erstellen
		onePath(findOnePath), triangle2(findTriangle2), square3(findSquare3),
		oneSquare(findOneSquare), triangleGE05(findTriangleGE05),
		//Die Heuristiken aus dem Paper erstellen
		heuristic1({&onePath, &triangle2, &square3}),
		heuristic2({&onePath, &oneSquare, &triangle2}),
		heuristic3({&onePath, &triangleGE05, &square3}),
		heuristic4({&onePath, &oneSquare, &triangleGE05}),
		tsp(tsp), lpData(lpData), tolerance(1e-5) {}


CutGenerator::CutStatus CombCutGen::validate(LinearProgram& lp, const std::vector<double>& solution,
											 CutStatus currentStatus) {
	//Dieser CutGen ist langsam, deshalb wird er nur genutzt, wenn die anderen Generatoren keine Cuts erzeugt haben
	if (currentStatus != CutGenerator::valid) {
		return CutGenerator::valid;
	}

	//Verletzte Constraints nach "Verletztheit" sortiert speichern
	std::set<tsp_util::ConstraintWithSlack> allConstrs;
	for (const CombHeuristic& ch:{heuristic1, heuristic2, heuristic3, heuristic4}) {
		std::vector<CombHeuristic::Comb> combs = ch.findViolatedCombs(lpData, solution);
		for (CombHeuristic::Comb& c:combs) {
			c.simplify(lpData, solution);
			LinearProgram::Constraint constr = getContraintFor(c);
			double lhs = constr.evalLHS(solution);
			double rhs = constr.getRHS();
			//Kann wegen Rundung auftreten
			if (tolerance.less(rhs, lhs)) {
				allConstrs.insert({constr, rhs - lhs});
			} else if (lemon::Tolerance<double>(0.5).less(lhs, rhs)) {
				/*
				 * Wenn die Ungleichung aber "sehr erfüllt" ist (viel slack), liegt dies wahrscheinlich an einem Fehler
				 * im Algorithmus
				 */
				std::cerr << "Found non-violated combs constraint: " << lhs << " vs " << rhs << std::endl;
			}
		}
	}

	if (!allConstrs.empty()) {
		size_t sumNZ = 0;
		//Constraint nach absteigender Verletztheit hinzufügen, bis zu viele Nonzeroes zum LP hinzugefügt wurden
		auto firstAboveMax = allConstrs.begin();
		while (firstAboveMax != allConstrs.end() &&
			   static_cast<city_id>(sumNZ) < tsp.getCityCount() * tsp.getCityCount()) {
			sumNZ += firstAboveMax->constraint.getNonzeroes().size();
			++firstAboveMax;
		}
		lp.addConstraints(allConstrs.begin(), firstAboveMax);
		return CutGenerator::maybe_recalc;
	} else {
		return CutGenerator::valid;
	}
}

//Gibt die entsprechende Kamm-Constraint zurück
LinearProgram::Constraint CombCutGen::getContraintFor(const CombHeuristic::Comb& c) {
	std::vector<double> coeffs(lpData.getVariableCount(), 0);
	std::vector<variable_id> indices;
	size_t rhs = lpData.sparserInducedSum(c.handle, coeffs, indices) + 1;
	for (const std::vector<city_id>& tooth:c.teeth) {
		rhs += lpData.sparserInducedSum(tooth, coeffs, indices);
	}
	rhs -= (c.teeth.size() + 1) / 2;
	return LinearProgram::Constraint::fromDense(indices, coeffs, LinearProgram::less_eq, rhs);
}

//Transformation 1 aus dem Paper
ContractionRule::Contraction CombCutGen::findOnePath(const Graph& g, const std::vector<Graph::Node>& possibleNodes,
													 const std::vector<Graph::Edge>&,
													 const Graph::EdgeMap<double>& costs,
													 const Graph::NodeMap<bool>& used) {
	lemon::Tolerance<double> tolerance;
	for (Graph::Node middle : possibleNodes) {
		size_t degree = 0;
		std::array<Graph::Node, 2> neighbors;
		for (Graph::OutArcIt it(g, middle); it != lemon::INVALID; ++it) {
			++degree;
			Graph::Node target = g.target(it);
			if (degree > 2 || tolerance.different(costs[it], 1) || used[target]) {
				//Mehr als 2 Kanten, eine Kante mit Wert !=1 oder schon genutzt->Kann nicht als mittlerer Knoten genutzt werden
				degree = 10;
				break;
			}
			neighbors[degree - 1] = target;
		}
		if (degree == 2) {
			assert(neighbors[0] != neighbors[1]);
			return {{middle, neighbors[0]},
					{neighbors[1]}};
		}
	}
	return {};
}

//Transformation 2 aus dem Paper
ContractionRule::Contraction CombCutGen::findTriangle2(const Graph& g, const std::vector<Graph::Node>&,
													   const std::vector<Graph::Edge>& possibleOneEdges,
													   const Graph::EdgeMap<double>& costs,
													   const Graph::NodeMap<bool>&) {
	lemon::Tolerance<double> tolerance;
	//Über alle möglichen Kanten {u, v} iterieren
	for (Graph::Edge oneEdge:possibleOneEdges) {
		//Kosten der an u anliegenden Kante zum entsprechenden Knoten
		Graph::NodeMap<double> adjCosts(g, 0);
		for (Graph::OutArcIt it(g, g.u(oneEdge)); it != lemon::INVALID; ++it) {
			adjCosts[g.target(it)] = costs[it];
		}
		//Über alle möglichen Knoten w iterieren
		for (Graph::OutArcIt it(g, g.v(oneEdge)); it != lemon::INVALID; ++it) {
			if (oneEdge != it && !tolerance.different(adjCosts[g.target(it)] + costs[it], 1)) {
				assert(g.target(it) != g.v(oneEdge));
				return {{g.u(oneEdge), g.v(oneEdge)},
						{g.target(it)}};
			}
		}
	}
	return {};
}

ContractionRule::Contraction CombCutGen::findSquare3(const Graph& g, const std::vector<Graph::Node>& possibleNodes,
													 const std::vector<Graph::Edge>&,
													 const Graph::EdgeMap<double>& costs,
													 const Graph::NodeMap<bool>& used) {
	lemon::Tolerance<double> tolerance;
	for (Graph::Node first:possibleNodes) {
		//Kosten der Kanten an first
		Graph::NodeMap<double> adjCosts1(g, 0);
		for (Graph::OutArcIt it(g, first); it != lemon::INVALID; ++it) {
			adjCosts1[g.target(it)] = costs[it];
		}

		//Mindestens ein Knoten im "Quadrat" muss mit 2 oder mehr der anderen Knoten im Quadrat benachbart sein
		for (Graph::OutArcIt itSecond(g, first); itSecond != lemon::INVALID; ++itSecond) {
			Graph::Node second = g.target(itSecond);
			if (used[second]) {
				continue;
			}
			//Kosten der Kanten an second
			Graph::NodeMap<double> adjCosts2(g, 0);
			for (Graph::OutArcIt it(g, second); it != lemon::INVALID; ++it) {
				adjCosts2[g.target(it)] = costs[it];
			}

			//Alle weiteren Nachbarn von 1 kommen als dritter Knoten in Frage
			Graph::OutArcIt itThird = itSecond;
			++itThird;
			for (; itThird != lemon::INVALID; ++itThird) {
				Graph::Node third = g.target(itThird);
				if (used[third]) {
					continue;
				}
				//Kosten der Kanten an third
				Graph::NodeMap<double> adjCosts3(g, 0);
				for (Graph::OutArcIt it(g, third); it != lemon::INVALID; ++it) {
					adjCosts3[g.target(it)] = costs[it];
				}

				for (Graph::Node fourth:possibleNodes) {
					if (fourth == first || fourth == second || fourth == third) {
						continue;
					}
					double indSum = adjCosts1[second] + adjCosts1[third] + adjCosts1[fourth]
									+ adjCosts2[third] + adjCosts2[fourth]
									+ adjCosts3[fourth];
					if (!tolerance.different(indSum, 3)) {
						return {{first, second, third, fourth}};
					}
				}
			}
		}
	}
	return {};
}

//Transformation 4
ContractionRule::Contraction CombCutGen::findOneSquare(const Graph& g, const std::vector<Graph::Node>&,
													   const std::vector<Graph::Edge>& possibleOneEdges,
													   const Graph::EdgeMap<double>& costs,
													   const Graph::NodeMap<bool>&) {
	//Wir brauchen 2 Kanten mit Wert 1
	if (possibleOneEdges.size() < 2) {
		return {};
	}
	lemon::Tolerance<double> tolerance;
	//Über alle möglichen Kombinationen von 1-Kanten iterieren
	for (size_t i = 0; i < possibleOneEdges.size() - 1; ++i) {
		Graph::Edge uw = possibleOneEdges[i];
		Graph::Node u = g.u(uw);
		Graph::Node w = g.v(uw);
		for (size_t i2 = i + 1; i2 < possibleOneEdges.size(); ++i2) {
			Graph::Edge vx = possibleOneEdges[i2];
			Graph::Node v = g.u(vx);
			Graph::Node x = g.v(vx);
			if (v == u || v == w || x == u || x == w) {
				continue;
			}

			//Erste verbindende Kante finden
			for (Graph::OutArcIt it(g, v); it != lemon::INVALID; ++it) {
				Graph::Node target = g.target(it);
				if (target == u || target == w) {
					//Ende der zweiten verbindenden Kante
					Graph::Node otherTarget = target == u ? w : u;
					for (Graph::OutArcIt it2(g, x); it2 != lemon::INVALID; ++it2) {
						if (g.target(it2) == otherTarget) {
							if (!tolerance.different(costs[it] + costs[it2], 1)) {
								return {{u, w},
										{v, x}
								};
							} else {
								break;
							}
						}
					}
					break;
				}
			}
		}

	}
	return {};
}

//Transformation 5
ContractionRule::Contraction CombCutGen::findTriangleGE05(const Graph& g, const std::vector<Graph::Node>&,
														  const std::vector<Graph::Edge>& possibleOneEdges,
														  const Graph::EdgeMap<double>& costs,
														  const Graph::NodeMap<bool>&) {
	for (Graph::Edge oneEdge:possibleOneEdges) {
		Graph::Node u = g.u(oneEdge);
		Graph::Node v = g.v(oneEdge);
		//Kosten der Kanten an u
		Graph::NodeMap<double> adjUCosts(g, 0);
		for (Graph::OutArcIt it(g, u); it != lemon::INVALID; ++it) {
			adjUCosts[g.target(it)] = costs[it];
		}
		for (Graph::OutArcIt it(g, v); it != lemon::INVALID; ++it) {
			Graph::Node w = g.target(it);
			double uwCost = adjUCosts[w];
			if (oneEdge != it && uwCost > 0 && uwCost + costs[it] >= 0.5) {
				assert(w != v && w != u);
				return {{u, v},
						{w}};
			}
		}
	}
	return {};
}