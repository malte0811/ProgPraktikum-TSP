#include <two_matching_cut_gen.hpp>
#include <lemon/unionfind.h>
#include <union_find.hpp>
#include <cmath>
#include <cassert>

TwoMatchingCutGen::TwoMatchingCutGen(const TSPInstance& inst, bool contract)
		: tsp(inst), enableContraction(contract), tolerance(1e-5) {}

CutGenerator::CutStatus TwoMatchingCutGen::validate(LinearProgram& lp, const std::vector<double>& solution) {
	Graph workGraph;
	std::vector<Graph::Node> origToWork(static_cast<size_t>(tsp.getCityCount()));
	ContractionMap workToOrig(workGraph);
	for (city_id i = 0; i<tsp.getCityCount(); ++i) {
		Graph::Node newNode = workGraph.addNode();
		origToWork[i] = newNode;
		workToOrig[newNode] = {i};
	}
	Graph::NodeMap<bool> odd(workGraph);
	Graph::EdgeMap <variable_id> toVariable(workGraph);
	Graph::EdgeMap<double> c(workGraph);
	std::vector<variable_id> oneEdges;
	for (city_id lower = 0; lower<tsp.getCityCount()-1; ++lower) {
		for (city_id higher = lower+1; higher<tsp.getCityCount(); ++higher) {
			variable_id varId = tsp.getVariable(higher, lower);
			if (tolerance.positive(solution[varId])) {
				Graph::Node endU = origToWork[lower];
				Graph::Node endV = origToWork[higher];
				if (tolerance.less(solution[varId], 1)) {
					Graph::Edge eWork = workGraph.addEdge(endU, endV);
					toVariable[eWork] = varId;
					c[eWork] = solution[varId];
				} else {
					odd[endU] = !odd[endU];
					odd[endV] = !odd[endV];
					oneEdges.push_back(varId);
				}
			}
		}
	}
	if (enableContraction) {
		contractPaths(workGraph, odd, c, workToOrig);
	}
	//TODO explain why I don't have z (c is always 0, cDash always inf, so the edges are irrelevant)
	std::vector<Blossom> allMin;
	lemma1220(workGraph, allMin, odd, c);
	if (!allMin.empty()) {
		std::vector<variable_id> indices;
		std::vector<int> constrStarts;
		std::vector<double> rhs;
		for (Blossom& min:allMin) {
			std::vector<bool> isInX(tsp.getCityCount());
			size_t sizeXTrue = 0;
			for (Graph::Node contracted:min.x) {
				for (city_id orig:workToOrig[contracted]) {
					isInX[orig] = true;
				}
				sizeXTrue += workToOrig[contracted].size();
			}
			std::vector<bool> fTSP(tsp.getEdgeCount());
			std::vector<variable_id> prelimF;
			for (variable_id e:oneEdges) {
				city_id endU = tsp.getLowerEnd(e);
				city_id endV = tsp.getHigherEnd(e);
				if (isInX[endU]!=isInX[endV]) {
					fTSP[e] = true;
					prelimF.push_back(e);
				}
			}
			for (Graph::Edge e:min.f) {
				fTSP[toVariable[e]] = true;
				prelimF.push_back(toVariable[e]);
			}
			std::vector<variable_id> incidentF(tsp.getEdgeCount(), LinearProgram::invalid_variable);
			size_t sizeF = prelimF.size();
			for (variable_id e:prelimF) {
				city_id ends[2] = {tsp.getLowerEnd(e), tsp.getHigherEnd(e)};
				for (city_id end:ends) {
					if (incidentF[end]==LinearProgram::invalid_variable) {
						incidentF[end] = e;
					} else {
						variable_id oldEdge = incidentF[end];
						incidentF[tsp.getLowerEnd(oldEdge)] = LinearProgram::invalid_variable;
						incidentF[tsp.getHigherEnd(oldEdge)] = LinearProgram::invalid_variable;
						incidentF[tsp.getLowerEnd(e)] = LinearProgram::invalid_variable;
						incidentF[tsp.getHigherEnd(e)] = LinearProgram::invalid_variable;
						fTSP[oldEdge] = false;
						fTSP[e] = false;
						isInX[end] = !isInX[end];
						if (isInX[end]) {
							++sizeXTrue;
						} else {
							--sizeXTrue;
						}
						sizeF -= 2;
						break;
					}
				}
			}
			/*
			 * Mit |F|==1 ist auch die Subtour-Constraint für X verletzt und impliziert die
			 * 2-Matching-Constraint
			 */
			constrStarts.push_back(static_cast<int>(indices.size()));
			if (sizeF>1) {
				for (variable_id e:prelimF) {
					if (fTSP[e]) {
						indices.push_back(e);
					}
				}
			}
			std::vector<city_id> xElements;
			const bool valForX = sizeXTrue<tsp.getCityCount()/2;
			for (city_id i = 0; i<tsp.getCityCount(); ++i) {
				if (isInX[i]==valForX) {
					for (city_id other:xElements) {
						indices.push_back(tsp.getVariable(i, other));
					}
					xElements.push_back(i);
				}
			}
			if (sizeF>1) {
				rhs.push_back(static_cast<size_t>(xElements.size()+sizeF/2));
			} else {
				rhs.push_back(xElements.size()-1);
			}
		}
		lp.addConstraints(indices, std::vector<double>(indices.size(), 1), rhs, constrStarts,
						  std::vector<LinearProgram::CompType>(rhs.size(), LinearProgram::less_eq));
		return CutGenerator::maybe_recalc;
	} else {
		return CutGenerator::valid;
	}
}

void TwoMatchingCutGen::lemma1220(const Graph& graph, std::vector<Blossom>& out, const Graph::NodeMap<bool>& odd,
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
	Graph::NodeMap <size_t> childCount(graph);
	size_t nodeCount = 0;
	for (Graph::NodeIt it(graph); it!=lemon::INVALID; ++it) {
		Graph::Node pred = gh.predNode(it);
		if (pred!=lemon::INVALID) {
			++childCount[pred];
		}
		++nodeCount;
	}
	Graph::NodeMap <size_t> ufMap(graph);
	UnionFind components(nodeCount);
	std::vector<Graph::Node> leaves;
	{
		size_t nextIndex = 0;
		for (Graph::NodeIt it(graph); it!=lemon::INVALID; ++it, ++nextIndex) {
			if (childCount[it]==0) {
				leaves.push_back(it);
			}
			ufMap[it] = nextIndex;
		}
	}
	while (!leaves.empty()) {
		Graph::Node currentNode = leaves.back();
		leaves.pop_back();
		Graph::Node pred = gh.predNode(currentNode);
		if (pred==lemon::INVALID) {
			continue;
		}
		const size_t xIndex = components.find(ufMap[currentNode]);
		double cutCost = gh.minCutValue(currentNode, pred);
		if (tolerance.less(cutCost, 1)) {
			Blossom curr;
			Graph::Edge minDiffEdge;
			double minDiffVal = std::numeric_limits<double>::max();
			bool hasMinEdge = false;
			size_t minIndex = std::numeric_limits<size_t>::max();
			for (Graph::EdgeIt eIt(graph); eIt!=lemon::INVALID; ++eIt) {
				bool uInX = components.find(ufMap[graph.u(eIt)])==xIndex;
				bool vInX = components.find(ufMap[graph.v(eIt)])==xIndex;
				if (uInX!=vInX) {
					bool inF = c[eIt]>0.5;
					if (inF) {
						curr.f.push_back(eIt);
					}
					if (std::abs(c[eIt]-0.5)<std::abs(minDiffVal)) {
						minDiffEdge = eIt;
						minDiffVal = c[eIt]-0.5;
						hasMinEdge = true;
						if (inF) {
							minIndex = curr.f.size()-1;
						} else {
							minIndex = std::numeric_limits<size_t>::max();
						}
					}
				}
			}
			unsigned xAndTDash = 0;
			for (Graph::NodeIt nIt(graph); nIt!=lemon::INVALID; ++nIt) {
				if (components.find(ufMap[nIt])==xIndex) {
					if ((adjacentEDash[nIt]%2==1)!=odd[nIt]) {
						++xAndTDash;
					}
					curr.x.push_back(nIt);
				}
			}
			bool valid = true;
			if (xAndTDash%2==0) {
				if (!hasMinEdge) {
					valid = false;
				} else if (minIndex<curr.f.size()) {
					cutCost += 2*minDiffVal;
					curr.f[minIndex] = curr.f.back();
					curr.f.pop_back();
				} else {
					cutCost -= 2*minDiffVal;
					curr.f.push_back(minDiffEdge);
				}
			}
			if (valid && tolerance.less(cutCost, 1)) {
				out.push_back(curr);
			}
		}
		components.mergeRoots(xIndex, components.find(ufMap[pred]));
		--childCount[pred];
		if (childCount[pred]==0) {
			leaves.push_back(pred);
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
	std::vector<city_id>& contractedSet = toOrig[remainingNode];
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