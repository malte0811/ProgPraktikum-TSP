#include <two_matching_cut_gen.hpp>

TwoMatchingCutGen::TwoMatchingCutGen(const TSPInstance& inst, bool contract)
		: tsp(inst), enableContraction(contract) {}

bool TwoMatchingCutGen::validate(LinearProgram& lp, const std::vector<double>& solution) {
	Graph workGraph;
	Graph::NodeMap <Graph::Node> origToWork(tsp.getGraph());
	Graph::NodeMap <std::vector<Graph::Node>> workToOrig(workGraph);
	for (Graph::NodeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		Graph::Node newNode = workGraph.addNode();
		origToWork[it] = newNode;
		workToOrig[newNode] = {it};
	}
	Graph::NodeMap<bool> odd(workGraph);
	Graph::EdgeMap <variable_id> vars(workGraph);
	Graph::EdgeMap<double> c(workGraph);
	std::vector<Graph::Edge> oneEdges;
	for (Graph::EdgeIt it(tsp.getGraph()); it!=lemon::INVALID; ++it) {
		variable_id varId = tsp.getVariable(it);
		if (tolerance.positive(solution[varId])) {
			Graph::Node endU = origToWork[tsp.getGraph().u(it)];
			Graph::Node endV = origToWork[tsp.getGraph().v(it)];
			if (tolerance.less(solution[varId], 1)) {
				Graph::Edge eWork = workGraph.addEdge(endU, endV);
				vars[eWork] = varId;
				c[eWork] = solution[varId];
			} else {
				odd[endU] = !odd[endU];
				odd[endV] = !odd[endV];
				oneEdges.push_back(it);
			}
		}
	}
	if (enableContraction) {
		contractPaths(workGraph, odd, c, vars, workToOrig);
	}
	//TODO explain why I don't have z (c is always 0, cDash always inf, so the edges are irrelevant)
	//TODO is that still true after contracting?
	std::vector<XandF> allMin;
	lemma1220(workGraph, allMin, odd, c);
	if (!allMin.empty() && tolerance.less(allMin.front().cost, 1)) {
		for (XandF& min:allMin) {
			std::vector<int> inSum;
			std::vector<city_id> xElements;
			std::vector<bool> isInX(tsp.getSize());
			for (Graph::Node contracted:min.x) {
				for (Graph::Node orig:workToOrig[contracted]) {
					variable_id uId = tsp.getCity(orig);
					isInX[uId] = true;
					for (int vId : xElements) {
						inSum.push_back(tsp.getVariable(uId, vId));
					}
					xElements.push_back(uId);
				}
			}
			size_t sizeF = min.f.size();
			for (Graph::Edge e:oneEdges) {
				Graph::Node endU = tsp.getGraph().u(e);
				Graph::Node endV = tsp.getGraph().v(e);
				city_id idU = tsp.getCity(endU);
				city_id idV = tsp.getCity(endV);
				if (isInX[idU]!=isInX[idV]) {
					inSum.push_back(tsp.getVariable(idU, idV));
					++sizeF;
				}
			}
			for (Graph::Edge e:min.f) {
				inSum.push_back(vars[e]);
			}
			double actualValue = 0;
			unsigned hash = 0;
			for (int a:inSum) hash = 31*hash+a, actualValue += solution[a];
			lp.addConstraint(inSum, std::vector<double>(inSum.size(), 1),
							 static_cast<int>(xElements.size()+sizeF/2),
							 LinearProgram::less_eq);
		}
		return false;
	} else {
		//if (enableContraction) {
		//	TwoMatchingCutGen verifier(tsp, false);
		//	if (!verifier.validate(lp, solution)) {
		//		std::cout << "ERROR, we found: " << allMin.size();
		//		if (!allMin.empty()) {
		//			std::cout << " with weight " << allMin[0].cost;
		//		}
		//		std::cout << std::endl;
		//	} else {
		//		//std::cout << "OK" << std::endl;
		//	}
		//}
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
									  Graph::EdgeMap <variable_id>& vars,
									  Graph::NodeMap <std::vector<Graph::Node>>& toOrig) {
	//std::vector<XandF> allMin;
	//lemma1220(g, allMin, odd, c);
	//double origVal = allMin.empty()?811:allMin.front().cost;
	bool contracted;
	do {
		contracted = false;
		for (Graph::NodeIt start(g); start!=lemon::INVALID; ++start) {
			if (odd[start]) {
				//TODO move to separate method
				Graph::NodeMap<bool> inPath(g);
				std::vector<Graph::Node> pathLeft, pathRight;
				double pathValue = -1;
				discoverPath(g, start, start, odd, inPath, pathLeft, c, toOrig, pathValue);
				if (pathLeft.empty()) {
					continue;
				}
				pathValue = 1-pathValue;
				discoverPath(g, start, pathLeft[0], odd, inPath, pathRight, c, toOrig, pathValue);
				if (pathLeft.size()+pathRight.size()<2) {
					continue;
				}
				//std::cout << "Contracting path of size " << pathLeft.size()+pathRight.size()+1 << std::endl << std::endl;
				std::vector<Graph::Node> toContract(pathLeft.begin(), pathLeft.end()-1);
				if (!pathRight.empty()) {
					toContract.insert(toContract.end(), pathRight.begin(), pathRight.end()-1);
				}
				Graph::Node remainingNode = pathLeft.back();
				if (remainingNode!=start) {
					toContract.push_back(start);
				}
				if (pathRight.back()!=pathLeft.back()) {
					inPath[pathRight.back()] = false;
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
					if (!g.valid(remove)) {
						std::cout << "Node is invalid!" << std::endl;
						continue;
					}
					if (remove==remainingNode) {
						std::cout << "Node is to remain!" << std::endl;
						continue;
					}
					if (odd[remove]) {
						resultOdd = !resultOdd;
					}
					contractedSet.insert(contractedSet.end(), toOrig[remove].begin(), toOrig[remove].end());
					g.erase(remove);
				}
				odd[remainingNode] = resultOdd;
				//std::vector<XandF> newAllMin;
				//lemma1220(g, newAllMin, odd, c);
				//double currVal = newAllMin.empty()?811:newAllMin.front().cost;
				//if (tolerance.less(origVal, 1) && tolerance.different(currVal, origVal)) {
				//	std::cout << "Last value: " << origVal << ", current: " << currVal << std::endl;
				//	origVal = currVal;
				//	allMin = newAllMin;
				//}
				contracted = true;
			}
		}
	} while (contracted);
}

void TwoMatchingCutGen::discoverPath(const Graph& graph, Graph::Node start, Graph::Node exclude,
									 const Graph::NodeMap<bool>& odd, Graph::NodeMap<bool>& visited,
									 std::vector<Graph::Node>& path, const Graph::EdgeMap<double>& c,
									 const Graph::NodeMap <std::vector<Graph::Node>>& toOrig, double& firstEdgeVal) {
	Graph::Node current = start;
	size_t degree;
	double edgeSum;
	double nextEdgeVal = firstEdgeVal;
	//std::cout << "Discovering path" << std::endl;
	do {
		if (current!=start) {
			path.push_back(current);
		}
		visited[current] = true;
		degree = 0;
		edgeSum = 0;
		//std::cout << "Next" << std::endl;
		Graph::Node nextNode;
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
			if (potential!=forbidden && (nextEdgeVal<0 || !tolerance.different(c[it], nextEdgeVal))) {
				nextNode = potential;
				if (tolerance.negative(nextEdgeVal)) {
					firstEdgeVal = c[it];
				}
				nextEdgeVal = c[it];
				//std::cout << c[it];
				//if (toOrig[current].size()>1) std::cout << " (Contracted)";
				//std::cout << std::endl;
			} else {
				//std::cout << c[it] << " (Forbidden)" << std::endl;
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
}
