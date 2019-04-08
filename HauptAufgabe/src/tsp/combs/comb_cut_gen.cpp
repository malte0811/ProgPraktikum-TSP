#include "comb_cut_gen.hpp"

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
			if (!tolerance.different(adjCosts[g.target(it)] + costs[it], 1)) {
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
	lemon::Tolerance<double> tolerance;
	for (size_t i = 0; i < possibleOneEdges.size() - 1; ++i) {
		Graph::Edge uw = possibleOneEdges[i];
		Graph::Node u = g.u(uw);
		Graph::Node w = g.v(uw);
		for (size_t i2 = i + 1; i2 < possibleOneEdges.size() - 1; ++i2) {
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
			if (!tolerance.less(adjCosts[g.target(it)] + costs[it], 0.5)) {
				return {{g.u(oneEdge), g.v(oneEdge)},
						{g.target(it)}};
			}
		}
	}
	return {};
}
