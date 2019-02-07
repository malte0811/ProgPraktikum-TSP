
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
	Graph::NodeMap<bool> odd(workGraph);
	Graph::EdgeMap<int> edgeToVar(workGraph);
	Graph::EdgeMap <variable_id> vars(workGraph);
	std::vector<Graph::Edge> oneEdges;
	for (Graph::EdgeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		variable_id varId = tsp.getVariable(it);
		if (solution[varId]>0) {
			Graph::Node endU = origToWork[tsp.getGraph().u(it)];
			Graph::Node endV = origToWork[tsp.getGraph().v(it)];
			if (solution[varId]==1) {
				odd[endU] = !odd[endU];
				odd[endV] = !odd[endV];
				oneEdges.push_back(it);
			} else {
				Graph::Edge eWork = workGraph.addEdge(endU, endV);
				vars[eWork] = varId;
				c[eWork] = solution[varId];
				cDash[eWork] = 1-solution[varId];
				edgeToVar[eWork] = varId;
			}
			totalVal[endU] += solution[varId];
			totalVal[endV] += solution[varId];
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
	lemma1220(allMin, odd, z);
	if (tolerance.less(allMin.front().cost, 1)) {
		for (XandF& min:allMin) {
			std::vector<int> inSum;
			std::vector<city_id> x;
			for (city_id uId = 0; uId<min.x.size(); ++uId) {
				if (min.x[uId]) {
					for (int vId : x) {
						inSum.push_back(tsp.getVariable(uId, vId));
					}
					x.push_back(uId);
				}
			}
			size_t sizeF = min.f.size();
			for (Graph::Edge e:oneEdges) {
				Graph::Node endU = tsp.getGraph().u(e);
				Graph::Node endV = tsp.getGraph().v(e);
				city_id idU = tsp.getCity(endU);
				city_id idV = tsp.getCity(endV);
				if (min.x[idU]!=min.x[idV]) {
					inSum.push_back(tsp.getVariable(idU, idV));
					++sizeF;
				}
			}
			for (Graph::Edge e:min.f) {
				inSum.push_back(vars[e]);
			}
			lp.addConstraint(inSum, std::vector<double>(inSum.size(), 1),
							 static_cast<int>(x.size()+sizeF/2),
							 LinearProgram::less_eq);
		}
		return false;
	} else {
		return true;
	}
}

void TwoMatchingCutGen::lemma1220(std::vector<XandF>& out, const Graph::NodeMap<bool>& odd, Graph::Node exclude) {
	Graph::EdgeMap<double> d(workGraph);
	Graph::NodeMap<int> adjacentEDash(workGraph);
	for (Graph::EdgeIt it(workGraph); it!=lemon::INVALID; ++it) {
		if (c[it]>cDash[it]) {
			++adjacentEDash[workGraph.v(it)];
			++adjacentEDash[workGraph.u(it)];
			d[it] = cDash[it];
		} else {
			d[it] = c[it];
		}
	}
	lemon::GomoryHu<Graph, Graph::EdgeMap<double>>
	gh(workGraph, d);
	gh.run();
	double minCost = std::numeric_limits<double>::max();
	Graph::NodeMap<bool> x(workGraph);
	Graph::EdgeMap<bool> f(workGraph);
	for (Graph::NodeIt it(workGraph); it!=lemon::INVALID; ++it) {
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
		bool hasMinEdge = false;
		for (Graph::EdgeIt eIt(workGraph); eIt!=lemon::INVALID; ++eIt) {
			if (x[workGraph.u(eIt)]!=x[workGraph.v(eIt)]) {
				f[eIt] = c[eIt]>cDash[eIt];
				if (!hasMinEdge || std::abs(c[eIt]-cDash[eIt])<std::abs(minDiffVal)) {
					minDiffEdge = eIt;
					minDiffVal = c[eIt]-cDash[eIt];
					hasMinEdge = true;
				}
			} else {
				f[eIt] = false;
			}
		}
		unsigned xAndTDash = 0;
		for (Graph::NodeIt nIt(workGraph); nIt!=lemon::INVALID; ++nIt) {
			if (x[nIt] && (adjacentEDash[nIt]%2==1)!=odd[nIt]) {
				++xAndTDash;
			}
		}
		if (xAndTDash%2==0) {
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
			XandF curr{std::vector<bool>(tsp.getSize()), std::vector<Graph::Edge>(0), cutCost};
			for (Graph::EdgeIt cp(workGraph); cp!=lemon::INVALID; ++cp) {
				if (f[cp]) {
					curr.f.push_back(cp);
				}
			}
			bool valForX = !x[exclude];
			for (Graph::NodeIt cp(workGraph); cp!=lemon::INVALID; ++cp) {
				if (x[cp]==valForX) {
					curr.x[tsp.getCity(workToOrig[cp])] = true;
				}
			}
			out.push_back(curr);
		}
	}
}