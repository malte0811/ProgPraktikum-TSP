
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
	Graph::EdgeMap <variable_id> vars(workGraph);
	for (Graph::EdgeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		variable_id varId = tsp.getVariable(it);
		if (solution[varId]>0) {
			Graph::Node uOrig = tsp.getGraph().u(it);
			Graph::Node vOrig = tsp.getGraph().v(it);
			Graph::Edge eWork = workGraph.addEdge(origToWork[uOrig], origToWork[vOrig]);
			vars[eWork] = varId;
			c[eWork] = solution[varId];
			cDash[eWork] = 1-solution[varId];
			edgeToVar[eWork] = varId;
			totalVal[origToWork[uOrig]] += solution[varId];
			totalVal[origToWork[vOrig]] += solution[varId];
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
	std::vector<XandF> allMin;
	lemma1220(workGraph, c, cDash, allMin, z);
	if (tolerance.less(allMin.front().cost, 1)) {
		for (XandF& min:allMin) {
			std::vector<int> inSum;
			for (size_t indexU = 1; indexU<min.x.size(); ++indexU) {
				for (size_t indexV = 0; indexV<indexU; ++indexV) {
					Graph::Node u = min.x[indexU];
					Graph::Node v = min.x[indexV];
					city_id uId = tsp.getCity(workToOrig[u]);
					city_id vId = tsp.getCity(workToOrig[v]);
					inSum.push_back(tsp.getVariable(uId, vId));
				}
			}
			for (Graph::Edge e:min.f) {
				inSum.push_back(vars[e]);
			}
			lp.addConstraint(inSum, std::vector<double>(inSum.size(), 1),
							 static_cast<int>(min.x.size()+min.f.size()/2),
							 LinearProgram::less_eq);
		}
		return false;
	} else {
		return true;
	}
}

void TwoMatchingCutGen::lemma1220(const Graph& graph, const Graph::EdgeMap<double>& c,
								  const Graph::EdgeMap<double>& cDash, std::vector<XandF>& out, Graph::Node exclude) {
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
	double minCost = std::numeric_limits<double>::max();
	Graph::NodeMap<bool> x(graph);
	Graph::EdgeMap<bool> f(graph);
	for (Graph::NodeIt it(graph); it!=lemon::INVALID; ++it) {
		Graph::Node pred = gh.predNode(it);
		if (pred==lemon::INVALID) {
			continue;
		}
		double cutCost = gh.minCutValue(it, pred);
		if (tolerance.less(minCost, cutCost)) {
			continue;
		}
		gh.minCutMap(it, pred, x);
		Graph::Edge minDiffEdge;
		double minDiffVal = std::numeric_limits<double>::max();
		for (Graph::EdgeIt eIt(graph); eIt!=lemon::INVALID; ++eIt) {
			if (x[graph.u(eIt)]!=x[graph.v(eIt)]) {
				f[eIt] = c[eIt]>cDash[eIt];
				if (std::abs(c[eIt]-cDash[eIt])<std::abs(minDiffVal)) {
					minDiffEdge = eIt;
					minDiffVal = c[eIt]-cDash[eIt];
				}
			} else {
				f[eIt] = false;
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
		if (!tolerance.less(minCost, cutCost)) {
			if (tolerance.less(cutCost, minCost)) {
				out.clear();
				minCost = cutCost;
			}
			XandF curr;
			curr.cost = cutCost;
			for (Graph::EdgeIt cp(graph); cp!=lemon::INVALID; ++cp) {
				if (f[cp]) {
					curr.f.push_back(cp);
				}
			}
			bool valForX = !x[exclude];
			for (Graph::NodeIt cp(graph); cp!=lemon::INVALID; ++cp) {
				if (x[cp]==valForX) {
					curr.x.push_back(cp);
				}
			}
			out.push_back(curr);
		}
	}
}