#include <subtour_cut_gen.hpp>
#include <lemon_fixes/nagamochi_ibaraki.h>
#include <cmath>
#include <lemon/dfs.h>

SubtourCutGen::SubtourCutGen(const TSPInstance& inst)
		: tsp(inst), origToWork(tsp.getGraph()), capacity(workGraph), minCut(workGraph, capacity) {
	for (Graph::NodeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		origToWork[it] = workGraph.addNode();
	}
	baseState.save(workGraph);
}

bool SubtourCutGen::validate(LinearProgram& lp, const std::vector<double>& solution) {
	baseState.restore();
	for (variable_id i = 0; i<solution.size(); ++i) {
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
	if (!tolerance.nonZero(capacity)) {
		Graph::NodeMap<bool> visited(workGraph);
		Graph::NodeIt it(workGraph);
		{
			//TODO comment (don't add constr for all comps)
			lemon::Dfs<Graph> dfs(workGraph);
			dfs.reachedMap(visited);
			dfs.run(it);
		}
		for (; it!=lemon::INVALID; ++it) {
			if (!visited[it]) {
				lemon::Dfs<Graph> dfs(workGraph);
				dfs.run(it);
				unsigned reached = 0;
				for (Graph::NodeIt it2(workGraph); it2!=lemon::INVALID; ++it2) {
					if (dfs.reached(it2)) {
						++reached;
						visited[it2] = true;
					}
				}
				std::vector<int> induced;
				for (Graph::EdgeIt eIt(tsp.getGraph()); eIt!=lemon::INVALID; ++eIt) {
					Graph::Node uOrig = tsp.getGraph().u(eIt);
					Graph::Node vOrig = tsp.getGraph().v(eIt);
					bool vInCut = dfs.reached(origToWork[vOrig]);
					bool uInCut = dfs.reached(origToWork[uOrig]);
					if (vInCut && uInCut) {
						induced.push_back(tsp.getVariable(eIt));
					}
				}
				lp.addConstraint(induced, std::vector<double>(induced.size(), 1), reached-1, LinearProgram::less_eq);
			}
		}
		return false;
	} else if (tolerance.less(capacity, 2)) {
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
		return false;
	} else {
		return true;
	}
}

