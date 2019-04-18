#include <comb_cut_gen.hpp>
#include <comb_heuristic.hpp>

CombCutGen::CombCutGen(const TSPInstance& tsp, const TspLpData& lpData) :
		onePath(contractOnePath), triangle2(contractTriangle2), square3(contractSquare3),
		oneSquare(contractOneSquare), triangleGE05(contractTriangleGE05),
		heuristic1({&onePath, &triangle2, &square3}),
		heuristic2({&onePath, &oneSquare, &triangle2}),
		heuristic3({&onePath, &triangleGE05, &square3}),
		heuristic4({&onePath, &oneSquare, &triangleGE05}),
		tsp(tsp), lpData(lpData) {}

CutGenerator::CutStatus CombCutGen::validate(LinearProgram& lp, const std::vector<double>& solution,
											 CutStatus currentStatus) {
	if (currentStatus != CutGenerator::valid) {
		return CutGenerator::valid;
	}
	std::vector<CombHeuristic::Comb> allCombs;
	for (const CombHeuristic& ch:{heuristic1, heuristic2, heuristic3, heuristic4}) {
		std::vector<CombHeuristic::Comb> combs = ch.findViolatedCombs(lpData, tsp, solution);
		for (CombHeuristic::Comb& c:combs) {
			c.simplify(lpData, solution);
		}
		allCombs.insert(allCombs.end(), combs.begin(), combs.end());
	}
	std::set<tsp_util::ConstraintWithSlack, tsp_util::CompareOrderedConstraint> allConstrs;
	for (const CombHeuristic::Comb& c:allCombs) {
		LinearProgram::Constraint constr = getContraintFor(c);
		//Kann wegen Rundung auftreten
		if (constr.isViolated(solution, lemon::Tolerance<double>(1e-4))) {
			allConstrs.insert({constr, constr.getRHS() - constr.evalLHS(solution)});
		}
	}
	if (!allConstrs.empty()) {
		std::vector<LinearProgram::Constraint> toAdd;
		size_t sumNZ = 0;
		for (const tsp_util::ConstraintWithSlack& pair:allConstrs) {
			sumNZ += pair.constraint.getNonzeroes().size();
			toAdd.push_back(pair.constraint);
			if (sumNZ > tsp.getCityCount() * tsp.getCityCount()) {
				break;
			}
		}
		std::cout << "Adding comb constraints with a total of " << sumNZ << " nonzeroes" << std::endl;
		lp.addConstraints(toAdd);
		return CutGenerator::maybe_recalc;
	} else {
		return CutGenerator::valid;
	}
}

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

ContractionRule::Contraction CombCutGen::contractOnePath(const Graph& g, const std::vector<Graph::Node>& possibleNodes,
														 const std::vector<Graph::Edge>& possibleOneEdges,
														 const Graph::EdgeMap<double>& costs,
														 const Graph::NodeMap<bool>& used) {
	lemon::Tolerance<double> tolerance;
	for (Graph::Node middle : possibleNodes) {
		size_t adjEdges = 0;
		std::array<Graph::Node, 2> neighbors;
		for (Graph::OutArcIt it(g, middle); it != lemon::INVALID; ++it) {
			++adjEdges;
			if (adjEdges > 2 || tolerance.different(costs[it], 1)) {
				adjEdges = 10;
				break;
			}
			Graph::Node target = g.target(it);
			if (used[target]) {
				adjEdges = 10;
				break;
			}
			neighbors[adjEdges - 1] = target;
		}
		if (adjEdges == 2) {
			assert(neighbors[0] != neighbors[1]);
			return {{middle, neighbors[0]},
					{neighbors[1]}};
		}
	}
	return {};
}

ContractionRule::Contraction CombCutGen::contractTriangle2(const Graph& g,
														   const std::vector<Graph::Node>& possibleNodes,
														   const std::vector<Graph::Edge>& possibleOneEdges,
														   const Graph::EdgeMap<double>& costs,
														   const Graph::NodeMap<bool>& used) {
	lemon::Tolerance<double> tolerance;
	for (Graph::Edge oneEdge:possibleOneEdges) {
		Graph::NodeMap<double> adjCosts(g, 0);
		for (Graph::OutArcIt it(g, g.u(oneEdge)); it != lemon::INVALID; ++it) {
			adjCosts[g.target(it)] = costs[it];
		}
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

ContractionRule::Contraction CombCutGen::contractSquare3(const Graph& g, const std::vector<Graph::Node>& possibleNodes,
														 const std::vector<Graph::Edge>& possibleOneEdges,
														 const Graph::EdgeMap<double>& costs,
														 const Graph::NodeMap<bool>& used) {
	//TODO das muss auch schöner gehen!
	lemon::Tolerance<double> tolerance;
	for (Graph::Node first:possibleNodes) {
		Graph::NodeMap<double> adjCosts1(g, 0);
		for (Graph::OutArcIt it(g, first); it != lemon::INVALID; ++it) {
			adjCosts1[g.target(it)] = costs[it];
		}
		for (Graph::OutArcIt it1(g, first); it1 != lemon::INVALID; ++it1) {
			Graph::Node second = g.target(it1);
			Graph::NodeMap<double> adjCosts2(g, 0);
			for (Graph::OutArcIt it(g, second); it != lemon::INVALID; ++it) {
				adjCosts2[g.target(it)] = costs[it];
			}
			if (used[second]) {
				continue;
			}
			Graph::OutArcIt it2 = it1;
			++it2;
			for (; it2 != lemon::INVALID; ++it2) {
				Graph::Node third = g.target(it2);
				Graph::NodeMap<double> adjCosts3(g, 0);
				for (Graph::OutArcIt it(g, third); it != lemon::INVALID; ++it) {
					adjCosts3[g.target(it)] = costs[it];
				}
				if (used[third]) {
					continue;
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

ContractionRule::Contraction CombCutGen::contractOneSquare(const Graph& g,
														   const std::vector<Graph::Node>& possibleNodes,
														   const std::vector<Graph::Edge>& possibleOneEdges,
														   const Graph::EdgeMap<double>& costs,
														   const Graph::NodeMap<bool>& used) {
	if (possibleOneEdges.size() < 2) {
		return {};
	}
	lemon::Tolerance<double> tolerance;
	for (size_t i = 0; i < possibleOneEdges.size() - 1; ++i) {
		Graph::Edge uw = possibleOneEdges[i];
		Graph::Node u = g.u(uw);
		Graph::Node w = g.v(uw);
		for (size_t i2 = i + 1; i2 < possibleOneEdges.size(); ++i2) {
			Graph::Edge vx = possibleOneEdges[i2];
			Graph::Node v = g.u(vx);
			Graph::Node x = g.v(vx);
			for (Graph::OutArcIt it(g, v); it != lemon::INVALID; ++it) {
				Graph::Node target = g.target(it);
				if (target != u && target != w) {
					continue;
				}
				Graph::Node otherTarget = target == u ? w : u;
				for (Graph::OutArcIt it2(g, x); it2 != lemon::INVALID; ++it2) {
					if (g.target(it2) == otherTarget && !tolerance.different(costs[it] + costs[it2], 1)) {
						return {{u, w},
								{v, x}};
					}
				}
			}
		}

	}
	return {};
}

ContractionRule::Contraction CombCutGen::contractTriangleGE05(const Graph& g,
															  const std::vector<Graph::Node>& possibleNodes,
															  const std::vector<Graph::Edge>& possibleOneEdges,
															  const Graph::EdgeMap<double>& costs,
															  const Graph::NodeMap<bool>& used) {
	lemon::Tolerance<double> tolerance;
	for (Graph::Edge oneEdge:possibleOneEdges) {
		Graph::NodeMap<double> adjCosts(g, 0);
		for (Graph::OutArcIt it(g, g.u(oneEdge)); it != lemon::INVALID; ++it) {
			adjCosts[g.target(it)] = costs[it];
		}
		for (Graph::OutArcIt it(g, g.v(oneEdge)); it != lemon::INVALID; ++it) {
			if (oneEdge != it && !tolerance.less(adjCosts[g.target(it)] + costs[it], 0.5)) {
				assert(g.target(it) != g.v(oneEdge));
				return {{g.u(oneEdge), g.v(oneEdge)},
						{g.target(it)}};
			}
		}
	}
	return {};
}