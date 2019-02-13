#include <two_matching_cut_gen.hpp>

TwoMatchingCutGen::TwoMatchingCutGen(const TSPInstance& inst, bool contract)
		: tsp(inst), enableContraction(contract), tolerance(1e-5) {}

bool TwoMatchingCutGen::validate(LinearProgram& lp, const std::vector<double>& solution) {
	Graph workGraph;
	Graph::NodeMap <Graph::Node> origToWork(tsp.getGraph());
	ContractionMap workToOrig(workGraph);
	for (Graph::NodeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		Graph::Node newNode = workGraph.addNode();
		origToWork[it] = newNode;
		workToOrig[newNode] = {it};
	}
	Graph::NodeMap<bool> odd(workGraph);
	Graph::EdgeMap <Graph::Edge> toOrigEdge(workGraph);
	Graph::EdgeMap<double> c(workGraph);
	std::vector<Graph::Edge> oneEdges;
	for (Graph::EdgeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		variable_id varId = tsp.getVariable(it);
		if (tolerance.positive(solution[varId])) {
			Graph::Node endU = origToWork[tsp.getGraph().u(it)];
			Graph::Node endV = origToWork[tsp.getGraph().v(it)];
			if (tolerance.less(solution[varId], 1)) {
				Graph::Edge eWork = workGraph.addEdge(endU, endV);
				toOrigEdge[eWork] = it;
				c[eWork] = solution[varId];
			} else {
				odd[endU] = !odd[endU];
				odd[endV] = !odd[endV];
				oneEdges.push_back(it);
			}
		}
	}
	if (enableContraction) {
		contractPaths(workGraph, odd, c, workToOrig);
	}
	//TODO explain why I don't have z (c is always 0, cDash always inf, so the edges are irrelevant)
	std::vector<XandF> allMin;
	lemma1220(workGraph, allMin, odd, c);
	if (!allMin.empty() && tolerance.less(allMin.front().cost, 1)) {
		for (XandF& min:allMin) {
			Graph::NodeMap<bool> isInX(tsp.getGraph());
			for (Graph::Node contracted:min.x) {
				for (Graph::Node orig:workToOrig[contracted]) {
					isInX[orig] = true;
				}
			}
			Graph::EdgeMap<bool> fTSP(tsp.getGraph());
			std::vector<Graph::Edge> prelimF;
			for (Graph::Edge e:oneEdges) {
				Graph::Node endU = tsp.getGraph().u(e);
				Graph::Node endV = tsp.getGraph().v(e);
				if (isInX[endU]!=isInX[endV]) {
					fTSP[e] = true;
					prelimF.push_back(e);
				}
			}
			for (Graph::Edge e:min.f) {
				fTSP[toOrigEdge[e]] = true;
				prelimF.push_back(toOrigEdge[e]);
			}
			Graph::NodeMap <Graph::Edge> incidentF(tsp.getGraph());
			for (Graph::NodeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
				incidentF[it] = lemon::INVALID;//TODO is this necessary?
			}
			for (Graph::Edge e:prelimF) {
				Graph::Node ends[2] = {tsp.getGraph().u(e), tsp.getGraph().v(e)};
				for (Graph::Node end:ends) {
					if (incidentF[end]==lemon::INVALID) {
						incidentF[end] = e;
					} else {
						Graph::Edge oldEdge = incidentF[end];
						fTSP[oldEdge] = false;
						incidentF[tsp.getGraph().u(oldEdge)] = lemon::INVALID;
						incidentF[tsp.getGraph().v(oldEdge)] = lemon::INVALID;
						incidentF[tsp.getGraph().u(e)] = lemon::INVALID;
						incidentF[tsp.getGraph().v(e)] = lemon::INVALID;
						fTSP[e] = false;
						isInX[end] = !isInX[end];
						break;
					}
				}
			}
			std::vector<variable_id> inSum;
			for (Graph::Edge e:prelimF) {
				if (fTSP[e]) {
					inSum.push_back(tsp.getVariable(e));
				}
			}
			const size_t sizeF = inSum.size();
			std::vector<city_id> xElements;
			for (Graph::NodeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
				if (isInX[it]) {
					city_id curr = tsp.getCity(it);
					for (city_id other:xElements) {
						inSum.push_back(tsp.getVariable(curr, other));
					}
					xElements.push_back(curr);
				}
			}
			size_t rhs = xElements.size()+sizeF/2;
			/*double actualValue = 0;
			unsigned hash = 0;
			for (int a:inSum) hash = 31*hash+a, actualValue += solution[a];
			if (actualValue<=rhs) {
				std::cout << "2-Matching cut-generator produced invalid (non-separating) cut! LHS " << actualValue << ", RHS " << rhs << std::endl;
				for (size_t id = 0;id<solution.size();++id) {
					if (tolerance.nonZero(solution[id])) {
						std::cout << id << ": " << solution[id] << "\n";
					}
				}
				std::cout << std::flush;
				return true;
			}*/
			lp.addConstraint(inSum, std::vector<double>(inSum.size(), 1),
							 rhs, LinearProgram::less_eq);
		}
		return false;
	} else {
		return true;
	}
}

void TwoMatchingCutGen::lemma1220(const Graph& graph, std::vector<XandF>& out, const Graph::NodeMap<bool>& odd,
								  const Graph::EdgeMap<double>& c) {
	Graph::EdgeMap<double> d(graph);
	Graph::NodeMap<int> adjacentEDash(graph);
	for (Graph::EdgeIt it(graph); it!=lemon::INVALID; ++it) {
		if (c[it]>0.5) {
			++adjacentEDash[graph.v(it)];
			++adjacentEDash[graph.u(it)];
			d[it] = 1-c[it];
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
		bool hasMinEdge = false;
		for (Graph::EdgeIt eIt(graph); eIt!=lemon::INVALID; ++eIt) {
			if (x[graph.u(eIt)]!=x[graph.v(eIt)]) {
				f[eIt] = c[eIt]>0.5;
				if (std::abs(c[eIt]-0.5)<std::abs(minDiffVal)) {
					minDiffEdge = eIt;
					minDiffVal = c[eIt]-0.5;
					hasMinEdge = true;
				}
			} else {
				f[eIt] = false;
			}
		}
		unsigned xAndTDash = 0;
		for (Graph::NodeIt nIt(graph); nIt!=lemon::INVALID; ++nIt) {
			if (x[nIt] && (adjacentEDash[nIt]%2==1)!=odd[nIt]) {
				++xAndTDash;
			}
		}
		if (xAndTDash%2==0) {
			if (!hasMinEdge) {
				continue;
			}
			if (f[minDiffEdge]) {
				cutCost += 2*minDiffVal;
				f[minDiffEdge] = false;
			} else {
				cutCost -= 2*minDiffVal;
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
			for (Graph::NodeIt cp(graph); cp!=lemon::INVALID; ++cp) {
				if (x[cp]) {
					curr.x.push_back(cp);
				}
			}
			out.push_back(curr);
		}
	}
}

void TwoMatchingCutGen::contractPaths(Graph& g, Graph::NodeMap<bool>& odd, Graph::EdgeMap<double>& c,
									  ContractionMap& toOrig) {
	bool contracted;
	do {
		contracted = false;
		for (Graph::NodeIt start(g); start!=lemon::INVALID && !contracted; ++start) {
			if (odd[start]) {
				contracted = findAndContractPath(g, start, toOrig, odd, c);
			}
		}
	} while (contracted);
}

bool TwoMatchingCutGen::findAndContractPath(Graph& g, Graph::Node start, ContractionMap& toOrig,
											Graph::NodeMap<bool>& odd,
											const Graph::EdgeMap<double>& c) {
	Graph::NodeMap<bool> inPath(g);
	double pathValue = -1;
	std::vector<Graph::Node> left = discoverPath(g, start, start, odd, inPath, c, pathValue);
	if (left.empty()) {
		return false;
	}
	pathValue = 1-pathValue;
	std::vector<Graph::Node> right = discoverPath(g, start, left[0], odd, inPath, c, pathValue);
	//TODO when can right be empty?
	if (left.size()+right.size()<2 || right.empty()) {
		return false;
	}
	std::vector<Graph::Node> toContract(left.begin(), left.end()-1);
	if (!right.empty()) {
		toContract.insert(toContract.end(), right.begin(), right.end()-1);
	}
	Graph::Node remainingNode = left.back();
	if (remainingNode!=start) {
		toContract.push_back(start);
	}
	if (right.back()!=left.back()) {
		inPath[right.back()] = false;
		for (Graph::Node pathNode:toContract) {
			std::vector<Graph::Edge> toMove;
			for (Graph::OutArcIt out(g, pathNode); out!=lemon::INVALID; ++out) {
				Graph::Node target = g.target(out);
				if (!inPath[target]) {
					toMove.push_back(out);
				}
			}
			for (Graph::Edge e:toMove) {
				if (g.u(e)==pathNode) {
					g.changeU(e, remainingNode);
				} else {
					g.changeV(e, remainingNode);
				}
			}
		}
	}
	bool resultOdd = odd[remainingNode];
	std::vector<Graph::Node>& contractedSet = toOrig[remainingNode];
	for (Graph::Node remove:toContract) {
		if (odd[remove]) {
			resultOdd = !resultOdd;
		}
		contractedSet.insert(contractedSet.end(), toOrig[remove].begin(), toOrig[remove].end());
		g.erase(remove);
	}
	odd[remainingNode] = resultOdd;
	return true;
}

std::vector<Graph::Node> TwoMatchingCutGen::discoverPath(const Graph& graph, Graph::Node start, Graph::Node exclude,
														 const Graph::NodeMap<bool>& odd, Graph::NodeMap<bool>& visited,
														 const Graph::EdgeMap<double>& c, double& firstEdgeVal) {
	std::vector<Graph::Node> path;
	Graph::Node current = start;
	size_t degree;
	double edgeSum;
	double nextEdgeVal = firstEdgeVal;
	do {
		if (current!=start) {
			path.push_back(current);
		}
		visited[current] = true;
		degree = 0;
		edgeSum = 0;
		Graph::Node nextNode = lemon::INVALID;
		Graph::Node forbidden;
		if (path.empty()) {
			forbidden = exclude;
		} else if (path.size()==1) {
			forbidden = start;
		} else {
			forbidden = path[path.size()-2];
		}
		for (Graph::OutArcIt it(graph, current); it!=lemon::INVALID; ++it) {
			Graph::Node potential = graph.target(it);
			if (potential!=forbidden && (tolerance.negative(nextEdgeVal) || !tolerance.different(c[it], nextEdgeVal))) {
				nextNode = potential;
				if (tolerance.negative(nextEdgeVal)) {
					firstEdgeVal = c[it];
				}
				nextEdgeVal = c[it];
			}
			++degree;
			edgeSum += c[it];
		}
		nextEdgeVal = 1-nextEdgeVal;
		current = nextNode;
	} while (graph.valid(current) && !visited[current] && odd[current] && degree==2 &&
			 !tolerance.different(edgeSum, 1));
	if (degree==2 && graph.valid(current) && !tolerance.different(edgeSum, 1)) {
		path.push_back(current);
		visited[current] = true;
	}
	return path;
}