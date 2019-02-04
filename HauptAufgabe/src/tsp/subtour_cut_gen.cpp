#include <subtour_cut_gen.hpp>
#include <lemon_fixes/nagamochi_ibaraki.h>
#include <cmath>

SubtourCutGen::SubtourCutGen(const TSPInstance& inst)
		: tsp(inst), origToWork(tsp.getGraph()), capacity(workGraph), minCut(workGraph, capacity) {
	for (Graph::NodeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		origToWork[it] = workGraph.addNode();
	}
	baseState.save(workGraph);
}

bool SubtourCutGen::validate(LinearProgram& lp, const std::vector<double>& solution) {
	baseState.restore();
	for (city_id i = 0; i<solution.size(); ++i) {
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
	if (tolerance.less(capacity, 2)) {
		Graph::NodeMap<bool> inCut(workGraph);
		minCut.minCutMap(inCut);
		city_id cutSize = 0;
		for (Graph::NodeIt it(workGraph); it!=lemon::INVALID; ++it) {
			if (inCut[it]) {
				++cutSize;
			}
		}
		std::vector<int> induced;
		//TODO figure out why this works so well
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

