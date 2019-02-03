
#include <two_matching_cut_gen.hpp>

TwoMatchingCutGen::TwoMatchingCutGen(const TSPInstance& inst)
		: tsp(inst), origToWork(inst.getGraph()), workToOrig(workGraph), c(workGraph), cDash(workGraph) {
	z = workGraph.addNode();
	for (Graph::NodeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		Graph::Node n = workGraph.addNode();
		origToWork[it] = n;
		workToOrig[n] = it;
		Graph::Edge e = workGraph.addEdge(z, origToWork[it]);
		cDash[e] = std::numeric_limits<double>::max();
	}
	basicState = Graph::Snapshot(workGraph);
}


bool TwoMatchingCutGen::validate(LinearProgram& lp, const std::vector<double>& solution) {
	basicState.restore();
	Graph::NodeMap<double> totalVal(workGraph);
	Graph::EdgeMap<int> edgeToVar(workGraph);
	Graph::EdgeMap <Graph::Edge> inWorkGraph(tsp.getGraph());
	for (Graph::EdgeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		variable_id varId = tsp.getVariable(it);
		if (solution[varId]>0) {
			Graph::Node uOrig = tsp.getGraph().u(it);
			Graph::Node vOrig = tsp.getGraph().v(it);
			Graph::Edge eWork = workGraph.addEdge(origToWork[uOrig], origToWork[vOrig]);
			inWorkGraph[it] = eWork;
			c[eWork] = solution[varId];
			cDash[eWork] = 1-solution[varId];
			edgeToVar[eWork] = varId;
			totalVal[workGraph.u(eWork)] += solution[varId];
			totalVal[workGraph.v(eWork)] += solution[varId];
		} else {
			inWorkGraph[it] = lemon::INVALID;
		}
	}
	for (Graph::IncEdgeIt it(workGraph, z); it!=lemon::INVALID; ++it) {
		//TODO is it always v/always u?
		Graph::Node otherEnd = workGraph.u(it);
		if (otherEnd==z) {
			otherEnd = workGraph.v(it);
		}
		c[it] = 2-totalVal[otherEnd];
	}
	XandF min(workGraph);
	lemma1220(workGraph, c, cDash, min);
	if (tolerance.less(min.cost, 1)) {
		std::vector<int> inSum;
		bool valForX = !min.x[z];
		if (!valForX) {
			//TODO better way of getting node count
			min.sizeX = workGraph.maxNodeId()+1-min.sizeX;
		}
		for (Graph::EdgeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
			Graph::Node uOrig = tsp.getGraph().u(it);
			Graph::Node vOrig = tsp.getGraph().v(it);
			if ((inWorkGraph[it]!=lemon::INVALID && min.f[inWorkGraph[it]]) ||
				(min.x[origToWork[uOrig]]==valForX && min.x[origToWork[vOrig]]==valForX)) {
				inSum.push_back(Graph::id(it));
			}
		}
		lp.addConstraint(inSum, std::vector<double>(inSum.size(), 1), min.sizeX+min.sizeF/2, LinearProgram::less_eq);
		return false;
	} else {
		return true;
	}
}

void TwoMatchingCutGen::lemma1220(const Graph& graph, const Graph::EdgeMap<double>& c,
								  const Graph::EdgeMap<double>& cDash, TwoMatchingCutGen::XandF& out) {
	Graph::EdgeMap<double> d(graph);
	Graph::NodeMap<int> adjacentEDash(graph);
	for (Graph::EdgeIt it(graph); it!=lemon::INVALID; ++it) {
		if (c[it]>cDash[it]) {
			++adjacentEDash[graph.v(it)];
			++adjacentEDash[graph.u(it)];
			d[it] = cDash[it];
		} else {
			d[it] = c[it];
		}
	}
	lemon::GomoryHu<Graph, Graph::EdgeMap<double>>
	gh(graph, d);
	gh.run();
	out.cost = std::numeric_limits<double>::max();
	for (Graph::NodeIt it(graph); it!=lemon::INVALID; ++it) {
		Graph::Node pred = gh.predNode(it);
		if (pred==lemon::INVALID) {
			continue;
		}
		Graph::NodeMap<bool> x(graph);
		gh.minCutMap(it, pred, x);
		Graph::EdgeMap<bool> f(graph);
		Graph::Edge minDiffEdge;
		double minDiffVal = std::numeric_limits<double>::max();
		double cutCost = 0;
		//TODO Geht cutCost<out.cost?
		for (Graph::EdgeIt eIt(graph); eIt!=lemon::INVALID; ++eIt) {
			if (x[graph.u(eIt)]!=x[graph.v(eIt)]) {
				if (c[eIt]>cDash[eIt]) {
					f[eIt] = true;
					cutCost += cDash[eIt];
				} else {
					cutCost += c[eIt];
				}
				if (std::abs(c[eIt]-cDash[eIt])<std::abs(minDiffVal)) {
					minDiffEdge = eIt;
					minDiffVal = c[eIt]-cDash[eIt];
				}
			}
		}
		unsigned xAndT = 0;
		for (Graph::NodeIt nIt(graph); nIt!=lemon::INVALID; ++nIt) {
			if (x[nIt] && adjacentEDash[nIt]%2==1) {
				++xAndT;
			}
		}
		if (xAndT%2==0) {
			if (f[minDiffEdge]) {
				cutCost += minDiffVal;
				f[minDiffEdge] = false;
			} else {
				cutCost -= minDiffVal;
				f[minDiffEdge] = true;
			}
		}
		if (cutCost<out.cost) {
			out.cost = cutCost;
			out.sizeF = out.sizeX = 0;
			for (Graph::EdgeIt cp(graph); cp!=lemon::INVALID; ++cp) {
				out.f[cp] = f[cp];
				out.sizeF += f[cp];
			}
			for (Graph::NodeIt cp(graph); cp!=lemon::INVALID; ++cp) {
				out.x[cp] = x[cp];
				out.sizeX += x[cp];
			}
		}
	}
}

TwoMatchingCutGen::XandF::XandF(const Graph& graph) : x(graph), f(graph) {}
