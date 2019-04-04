#include <one_tree.hpp>
#include <lemon/kruskal.h>
#include <numeric>

OneTree::OneTree(const TSPInstance& inst) : inst(&inst), potential(inst.getCityCount()-1, 0),
											g(inst.getCityCount() - 1), t(g, inTree), inTree(g, false), bestInTree(g),
											costs(g), tolerance(1e-3), lowerBound(0), lastV(inst.getCityCount() - 2) {
	for (FullGraph::EdgeIt it(g); it != lemon::INVALID; ++it) {
		sortedEdges.push_back(std::make_pair(it, costs[it]));
	}
}

void OneTree::run() {
	size_t iterations = 0;
	size_t sinceLastUpdate = 0;
	double currStepSize = 1;
	bool stop;
	const size_t maxIterations = 1000;
	const size_t maxNoUpdate = inst->getCityCount() / 10;
	const double factor = 0.99;
	do {
		generate1Tree();
		//Formel kommt von Martin Drees
		stop = updatePotential(currStepSize * 0.2 * costPotential / (double) g.nodeNum());
		++iterations;
		++sinceLastUpdate;
		if (lowerBound < costPotential) {
			lowerBound = costPotential;
			lemon::mapCopy(g, inTree, bestInTree);
			sinceLastUpdate = 0;
		}
		currStepSize *= factor;
	} while (!stop && sinceLastUpdate<maxNoUpdate && iterations < maxIterations);
	std::cout << iterations << ", " << sinceLastUpdate << ", " << currStepSize << std::endl;
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

void OneTree::generate1Tree() {
	//TODO This is a mess...
	//TODO wäre Prim besser? Gibt es in meiner Version von lemon leider nicht
	for (auto& pair:sortedEdges) {
		FullGraph::Edge e = pair.first;
		city_id endU = FullGraph::id(g.u(e)) + 1;
		city_id endV = FullGraph::id(g.v(e)) + 1;
		costs[e] = inst->getDistance(endU, endV) + potential[endU-1] + potential[endV-1];
		pair.second = costs[pair.first];
	}
	std::sort(sortedEdges.begin(), sortedEdges.end(), lemon::_kruskal_bits::PairComp<Sequence>());
	costPotential = lemon::_kruskal_bits::kruskal(g, sortedEdges, inTree);
	costOriginal = 0;
	double minEdgeCost = std::numeric_limits<double>::max();
	double secondMinEdgeCost = std::numeric_limits<double>::max();
	city_id minEdgeEnd = 0;
	city_id secondMinEdgeEnd = 0;
	double sumPotential = 0;
	for (city_id other = 1; other < inst->getCityCount(); ++other) {
		double eCost = inst->getDistance(0, other) + potential[other-1];
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
		sumPotential += potential[other-1];
	}
	oneNeighbors = {minEdgeEnd, secondMinEdgeEnd};
	costPotential += minEdgeCost + secondMinEdgeCost-2 * sumPotential;
	costOriginal += static_cast<cost_t>(minEdgeCost + secondMinEdgeCost
										- potential[minEdgeEnd-1] - potential[secondMinEdgeEnd-1]);
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

bool OneTree::updatePotential(double stepSize) {
	std::vector<int> v(inst->getCityCount()-2);
	for (int i = 1; i < inst->getCityCount() - 1; ++i) {
		FullGraph::Node curr = g(i - 1);
		int incidentCount = 0;
		for (Tree::IncEdgeIt it(t, curr); it != lemon::INVALID; ++it) {
			++incidentCount;
		}
		v[i-1] = incidentCount-2;
	}
	if (oneNeighbors.first < inst->getCityCount()-1) {
		v[oneNeighbors.first-1] += 1;
	}
	if (oneNeighbors.second < inst->getCityCount()-1) {
		v[oneNeighbors.second-1] += 1;
	}
	bool allZero = true;
	for (size_t i = 0;i<v.size();++i) {
		potential[i] += (0.6 * v[i] + 0.4 * lastV[i]) * stepSize;
		allZero &= v[i]==0;
		lastV[i] = v[i];
	}
	return allZero;
}

std::pair<city_id, city_id> OneTree::getOneNeighbors() const {
	return oneNeighbors;
}
