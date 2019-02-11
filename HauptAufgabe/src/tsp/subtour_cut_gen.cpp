#include <subtour_cut_gen.hpp>
#include <lemon_fixes/nagamochi_ibaraki.h>
#include <cmath>
#include <stack>

SubtourCutGen::SubtourCutGen(const TSPInstance& inst)
		: tsp(inst), origToWork(tsp.getGraph()), workToOrig(workGraph), capacity(workGraph),
		  minCut(workGraph, capacity) {
	for (Graph::NodeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		Graph::Node newNode = workGraph.addNode();
		origToWork[it] = newNode;
		workToOrig[newNode] = it;
	}
	baseState.save(workGraph);
}

bool SubtourCutGen::validate(LinearProgram& lp, const std::vector<double>& solution) {
	baseState.restore();
	baseState.save(workGraph);
	for (variable_id i = 0; i<static_cast<variable_id>(solution.size()); ++i) {
		if (tolerance.positive(solution[i])) {
			Graph::Edge inOrig = tsp.getEdge(i);
			Graph::Node uOrig = tsp.getGraph().u(inOrig);
			Graph::Node vOrig = tsp.getGraph().v(inOrig);
			Graph::Edge inWork = workGraph.addEdge(origToWork[uOrig], origToWork[vOrig]);
			capacity[inWork] = solution[i];
		}
	}
	minCut.run();
	double capacity = minCut.minCutValue();
	if (capacity==0) {
		addConnectivityConstraints(lp);
		return false;
	} else if (tolerance.less(capacity, 2)) {
		addCutConstraint(lp);
		return false;
	} else {
		return true;
	}
}

void SubtourCutGen::addConnectivityConstraints(LinearProgram& lp) {
	Graph::NodeMap<bool> visited(workGraph);
	bool firstComp = true;
	for (Graph::NodeIt it(workGraph); it!=lemon::INVALID; ++it) {
		if (!visited[it]) {
			std::stack<Graph::OutArcIt> stack;
			std::vector<city_id> currentComponent{tsp.getCity(workToOrig[it])};
			stack.push(Graph::OutArcIt(workGraph, it));
			visited[it] = true;
			while (!stack.empty()) {
				Graph::OutArcIt& current = stack.top();
				while (current!=lemon::INVALID) {
					Graph::Node neighbor = workGraph.target(current);
					if (!visited[neighbor]) {
						visited[neighbor] = true;
						stack.push(Graph::OutArcIt(workGraph, neighbor));
						currentComponent.push_back(tsp.getCity(workToOrig[neighbor]));
						break;
					}
					++current;
				}
				if (current==lemon::INVALID) {
					stack.pop();
				}
			}
			if (!firstComp) {
				std::vector<int> induced;
				for (size_t aId = 1; aId<currentComponent.size(); ++aId) {
					for (size_t bId = 0; bId<aId; ++bId) {
						induced.push_back(tsp.getVariable(currentComponent[aId], currentComponent[bId]));
					}
				}
				lp.addConstraint(induced, std::vector<double>(induced.size(), 1), currentComponent.size()-1,
								 LinearProgram::less_eq);
			} else {
				firstComp = false;
			}
		}
	}
}

void SubtourCutGen::addCutConstraint(LinearProgram& lp) {
	Graph::NodeMap<bool> inCut(workGraph);
	minCut.minCutMap(inCut);
	city_id cutSize = 0;
	for (Graph::NodeIt it(workGraph); it!=lemon::INVALID; ++it) {
		if (inCut[it]) {
			++cutSize;
		}
	}
	std::vector<int> induced;
	//TODO figure out why this works very well or not at all depending on edge ordering
	const bool cutVal = cutSize<tsp.getSize()/2;
	if (!cutVal) {
		cutSize = tsp.getSize()-cutSize;
	}
	for (Graph::EdgeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		Graph::Node uOrig = tsp.getGraph().u(it);
		Graph::Node vOrig = tsp.getGraph().v(it);
		bool vInCut = inCut[origToWork[vOrig]];
		bool uInCut = inCut[origToWork[uOrig]];
		if (vInCut==cutVal && uInCut==cutVal) {
			induced.push_back(tsp.getVariable(it));
		}
	}
	lp.addConstraint(induced, std::vector<double>(induced.size(), 1), cutSize-1, LinearProgram::less_eq);
}

