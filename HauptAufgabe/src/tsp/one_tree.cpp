#include <one_tree.hpp>
#include <lemon/kruskal.h>
#include <numeric>

OneTree::OneTree(const TSPInstance& inst) : inst(&inst), potential(inst.getCityCount(), 0),
											g(inst.getCityCount() - 1), t(g, inTree), inTree(g, false), bestInTree(g),
											costs(g), tolerance(1e-3) {}

void OneTree::run(CPXENVptr env) {
	LinearProgram lp(env, "1-Tree", LinearProgram::minimize);
	lp.addConstraint(LinearProgram::Constraint({}, {}, LinearProgram::equal, 1));
	lp.addConstraints(std::vector<LinearProgram::Constraint>(inst->getCityCount() - 2,
															 LinearProgram::Constraint({}, {}, LinearProgram::equal,
																					   0)));
	addInitialTrees(lp);
	bool stop;
	size_t iterations = 0;
	const size_t maxIterations = 10 * inst->getCityCount();
	do {
		LinearProgram::Solution sol = lp.solvePrimal();
		assert(sol.isValid());
		potential = sol.getShadowCosts();
		assert(!tolerance.different(potential[0], sol.getValue()));
		potential.push_back(0);//TODO das ist nur Faulheit
		stop = generate1Tree();
		if (!stop) {
			addVariable(lp);
		}
		++iterations;
	} while (!stop && iterations < maxIterations);
	lemon::mapCopy(g, bestInTree, inTree);
}

double OneTree::getCleanCost() const {
	return costOriginal;
}

const std::vector<double>& OneTree::getPotential() const {
	return potential;
}

const OneTree::Tree& OneTree::getTree() const {
	return t;
}

void OneTree::addInitialTrees(LinearProgram& lp) {
	for (city_id center = 1; center < inst->getCityCount(); ++center) {
		city_id deg2 = center % (inst->getCityCount() - 1) + 1;
		cost_t wheelCost = inst->getDistance(deg2, 0);
		for (city_id i = 0; i < inst->getCityCount(); ++i) {
			if (i != center) {
				wheelCost += inst->getDistance(center, i);
			}
		}
		std::vector<double> coeffs(inst->getCityCount() - 1, 1);
		std::vector<int> indices(inst->getCityCount() - 1);
		std::iota(std::begin(indices), std::end(indices), 0);
		coeffs[1] = 0;
		if (deg2 < inst->getCityCount() - 2) {
			coeffs[deg2 + 1] = 0;
		}
		if (center < inst->getCityCount() - 2) {
			coeffs[center + 1] = -(inst->getCityCount() - 3.0);
		}
		lp.addVariable(wheelCost, 0, 1, indices, coeffs);
	}
}

bool OneTree::generate1Tree() {
	for (FullGraph::EdgeIt it(g); it != lemon::INVALID; ++it) {
		city_id endU = FullGraph::id(g.u(it)) + 1;
		city_id endV = FullGraph::id(g.v(it)) + 1;
		costs[it] = inst->getDistance(endU, endV) + potential[endU] + potential[endV];
	}
	costPotential = lemon::kruskal(g, costs, inTree);
	costOriginal = 0;
	double minEdgeCost = std::numeric_limits<double>::max();
	double secondMinEdgeCost = std::numeric_limits<double>::max();
	city_id minEdgeEnd = 0;
	city_id secondMinEdgeEnd = 0;
	double sumPotential = 0;
	for (city_id other = 1; other < inst->getCityCount(); ++other) {
		double eCost = inst->getDistance(0, other) + potential[other];
		if (eCost <= minEdgeCost) {
			secondMinEdgeCost = minEdgeCost;
			secondMinEdgeEnd = minEdgeEnd;
			minEdgeCost = eCost;
			minEdgeEnd = other;
		} else if (eCost < secondMinEdgeCost) {
			secondMinEdgeCost = eCost;
			secondMinEdgeEnd = other;
		}
		Tree::Node curr = g(other - 1);
		for (Tree::OutArcIt it(t, curr); it != lemon::INVALID; ++it) {
			city_id targetId = t.id(t.target(it)) + 1;
			if (targetId > other) {
				costOriginal += inst->getDistance(other, targetId);
			}
		}
		sumPotential += potential[other];
	}
	oneNeighbors = {minEdgeEnd, secondMinEdgeEnd};
	costPotential += minEdgeCost + secondMinEdgeCost;
	costOriginal += static_cast<cost_t>(minEdgeCost + secondMinEdgeCost
										- potential[minEdgeEnd] - potential[secondMinEdgeEnd]);
	double currBound = costPotential - 2 * sumPotential;
	if (lowerBound < currBound) {
		lowerBound = currBound;
		lemon::mapCopy(g, inTree, bestInTree);
	}
	return !tolerance.less(currBound, potential[0]);
}

void OneTree::addVariable(LinearProgram& lp) {
	std::vector<double> coeffs(lp.getConstraintCount());
	std::vector<int> indices(lp.getConstraintCount());
	std::iota(std::begin(indices), std::end(indices), 0);
	coeffs[0] = 1;
	for (int i = 1; i < inst->getCityCount() - 1; ++i) {
		FullGraph::Node curr = g(i - 1);
		size_t incidentCount = 0;
		for (Tree::IncEdgeIt it(t, curr); it != lemon::INVALID; ++it) {
			++incidentCount;
		}
		coeffs[i] = 2 - static_cast<double>(incidentCount);
	}
	if (oneNeighbors.first < coeffs.size()) {
		coeffs[oneNeighbors.first] -= 1;
	}
	if (oneNeighbors.second < coeffs.size()) {
		coeffs[oneNeighbors.second] -= 1;
	}
	double sum = 0;
	for (size_t i = 0; i < coeffs.size(); ++i) {
		sum += coeffs[i] * potential[i];
	}
	lp.addVariable(costOriginal, 0, 1, indices, coeffs);
}

const OneTree::FullGraph& OneTree::getGraph() const {
	return g;
}

const lemon::FullGraph::EdgeMap<bool>& OneTree::getInTree() const {
	return inTree;
}

const lemon::GraphExtender<lemon::FullGraphBase>::EdgeMap<double>& OneTree::getReducedCosts() const {
	return costs;
}

double OneTree::getBound() {
	return lowerBound;
}
