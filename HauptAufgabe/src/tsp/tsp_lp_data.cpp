#include <tsp_lp_data.hpp>
#include <lemon/connectivity.h>

TspLpData::TspLpData(const TSPInstance& inst) : inst(inst), variableToEdge(inst.getEdgeCount()),
												edgeToVariable(inst.getCityCount() - 1),
												removalBound(inst.getEdgeCount(), -std::numeric_limits<double>::max()),
												shouldHaveRemoved(inst.getEdgeCount()) {
	variable_id currVar = 0;
	for (city_id higher = 1; higher < inst.getCityCount(); ++higher) {
		edgeToVariable[higher - 1].resize(higher);
		for (city_id lower = 0; lower < higher; ++lower) {
			variableToEdge[currVar] = {lower, higher};
			edgeToVariable[higher - 1][lower] = currVar;
			++currVar;
		}
	}
}

void TspLpData::setupLowerBounds() {
	OneTree oneTree(inst);
	oneTree.run();
	lemon::Bfs<OneTree::Tree> bfs(oneTree.getTree());
	bfs.run(oneTree.getGraph()(0));
	std::cout << "Lower bound is " << oneTree.getBound() << std::endl;
	assert(lemon::connected(oneTree.getTree()));
	const OneTree::FullGraph::EdgeMap<bool>& inTree = oneTree.getInTree();
	const OneTree::FullGraph::EdgeMap<double>& redCosts = oneTree.getReducedCosts();
	std::pair<city_id, city_id> oneNeighbors = oneTree.getOneNeighbors();
	OneTree::FullGraph::Node zero = oneTree.getGraph()(0);
	OneTree::FullGraph::Node neighbor1 = oneTree.getGraph()(oneNeighbors.first - 1);
	OneTree::FullGraph::Node neighbor2 = oneTree.getGraph()(oneNeighbors.second - 1);
	double expensiveOneEdge = std::max(redCosts[oneTree.getGraph().edge(zero, neighbor1)],
									   redCosts[oneTree.getGraph().edge(zero, neighbor2)]);
	for (variable_id i = 0; i < variableToEdge.size(); ++i) {
		Edge e = getEdge(i);
		if (e.first == 0 || e.second == 0) {
			//TODO removalBound[i] = oneTree.getBound()+1+redCost-expensiveOneEdge;
		} else {
			OneTree::FullGraph::Node endA = oneTree.getGraph()(e.first - 1);
			OneTree::FullGraph::Node endB = oneTree.getGraph()(e.second - 1);
			OneTree::FullGraph::Edge eGraph = oneTree.getGraph().edge(endA, endB);
			double redCost = redCosts[eGraph];
			if (inTree[eGraph])
				continue;
			double expensiveCycleEdge = -std::numeric_limits<double>::max();
			while (endA != endB) {
				OneTree::Tree::Arc next;
				if (bfs.dist(endA) > bfs.dist(endB)) {
					next = bfs.predArc(endA);
					endA = oneTree.getTree().source(next);
				} else {
					next = bfs.predArc(endB);
					endB = oneTree.getTree().source(next);
				}
				expensiveCycleEdge = std::max(expensiveCycleEdge, redCosts[next]);
			}
			removalBound[i] = oneTree.getBound() + 1 + redCost - expensiveCycleEdge;
		}

	}
}

std::vector<variable_id> TspLpData::removeVariables(coeff_t bound, const std::vector<value_t>& variables) {
	std::vector<bool> asBools(variables.size());
	for (size_t i = 0; i < asBools.size(); ++i) {
		asBools[i] = variables[i] != 0;
	}
	upperBound = TSPSolution(inst, asBools, *this);
	std::vector<variable_id> toRemove;
	variable_id newId = 0;
	for (variable_id oldId = 0; oldId < variableToEdge.size(); ++oldId) {
		//TODO tolernace?
		if (removalBound[oldId] > bound) {
			toRemove.push_back(oldId);
		} else {
			Edge e = variableToEdge[oldId];
			edgeToVariable[e.second - 1][e.first] = newId;
			++newId;
		}
	}
	for (auto it = toRemove.rbegin(); it != toRemove.rend(); ++it) {
		variable_id rem = *it;
		assert(removalBound[rem] > bound);
		removalBound.erase(removalBound.begin() + rem);
		Edge e = variableToEdge[rem];
		variableToEdge.erase(variableToEdge.begin() + rem);
		assert(edgeToVariable[e.second - 1][e.first] == rem);
		edgeToVariable[e.second - 1][e.first] = LinearProgram::invalid_variable;
	}
	for (variable_id i = 0; i < variableToEdge.size(); ++i) {
		Edge e = getEdge(i);
		assert(getVariable(e.first, e.second) == i);
	}
	return toRemove;
}


/**
 * Fügt Variablen für alle Kanten und Gradbedingungen für alle Knoten zum LP hinzu
 * @param lp das zu initialisierende LP
 */
void TspLpData::setupBasicLP(LinearProgram& lp) const {
	auto varCount = static_cast<size_t>(inst.getEdgeCount());
	std::vector<double> objCoeffs(varCount);
	for (variable_id i = 0; i < inst.getEdgeCount(); ++i) {
		objCoeffs[i] = getCost(i);
	}
	std::vector<double> lower(varCount, 0);
	std::vector<double> upper(varCount, 1);
	lp.addVariables(objCoeffs, lower, upper);
	std::vector<LinearProgram::Constraint> constrs;
	constrs.reserve(inst.getCityCount());
	for (city_id i = 0; i < inst.getCityCount(); ++i) {
		std::vector<variable_id> indices;
		for (city_id otherEnd = 0; otherEnd < inst.getCityCount(); ++otherEnd) {
			if (otherEnd != i)
				indices.push_back(getVariable(i, otherEnd));
		}
		constrs.emplace_back(indices, std::vector<double>(indices.size(), 1), LinearProgram::equal,
							 2);
	}
	lp.addConstraints(constrs);
}

const TSPSolution& TspLpData::getUpperBound() const {
	return upperBound;
}
