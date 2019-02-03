#include <subtour_cut_gen.hpp>
#include <lemon_fixes/nagamochi_ibaraki.h>
#include <cmath>

SubtourCutGen::SubtourCutGen(const lemon::SmartGraph& graph)
		: origGraph(graph), capacity(workGraph), minCut(workGraph, capacity) {}

bool SubtourCutGen::validate(LinearProgram& lp, const std::vector<double>& solution) {
	workGraph.clear();
	workGraph.reserveNode(origGraph.maxNodeId()+1);
	workGraph.reserveNode(origGraph.maxEdgeId()+1);
	for (node_id i = 0; i<=origGraph.maxNodeId(); ++i) {
		workGraph.addNode();
	}
	for (node_id i = 0; i<solution.size(); ++i) {
		if (tolerance.positive(solution[i])) {
			Graph::Edge inOrig = lemon::SmartGraph::edgeFromId(i);
			Graph::Edge inWork = workGraph.addEdge(origGraph.u(inOrig), origGraph.v(inOrig));
			capacity[inWork] = solution[i];
		}
	}
	minCut.run();
	double capacity = minCut.minCutValue();
	if (tolerance.less(capacity, 2)) {
		lemon::SmartGraph::NodeMap<bool> inCut(workGraph);
		minCut.minCutMap(inCut);
		std::vector<int> inducedEdges;
		inducedEdges.reserve(origGraph.maxEdgeId()+1);
		for (lemon::SmartGraph::EdgeIt it(origGraph); it!=lemon::INVALID; ++it) {
			if (inCut[origGraph.u(it)] && inCut[origGraph.v(it)]) {
				inducedEdges.push_back(lemon::SmartGraph::id(it));
			}
		}
		node_id cutSize = 0;
		for (lemon::SmartGraph::NodeIt it(origGraph); it!=lemon::INVALID; ++it) {
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

