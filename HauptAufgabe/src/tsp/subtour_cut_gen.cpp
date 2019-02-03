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
		std::vector<int> inducedEdges;
		for (Graph::EdgeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
			Graph::Node uOrig = tsp.getGraph().u(it);
			Graph::Node vOrig = tsp.getGraph().v(it);
			if (inCut[origToWork[uOrig]] && inCut[origToWork[vOrig]]) {
				inducedEdges.push_back(tsp.getVariable(it));
			}
		}
		city_id cutSize = 0;
		for (Graph::NodeIt it(workGraph); it!=lemon::INVALID; ++it) {
			if (inCut[it]) {
				++cutSize;
			}
		}
		lp.addConstraint(inducedEdges, std::vector<double>(inducedEdges.size(), 1), cutSize-1, LinearProgram::less_eq);
		return false;
	} else {
		return true;
	}
}

