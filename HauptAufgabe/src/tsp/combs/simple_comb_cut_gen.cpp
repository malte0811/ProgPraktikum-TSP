#include <simple_comb_cut_gen.hpp>
#include <tsp_utils.hpp>

SimpleCombCutGen::SimpleCombCutGen(const TSPInstance& inst, const TspLpData& lpData) : tsp(inst), tolerance(1e-5),
																					   lpData(lpData) {}

SimpleCombCutGen::BlockDecomposition SimpleCombCutGen::generateBlocks(const Graph& g) {
	Graph::EdgeMap<int> compMap(g);
	const int compCount = lemon::biNodeConnectedComponents(g, compMap);
	//std::cout << compCount << std::endl;
	Graph blockCutTree;
	BlockDecomposition::OriginalNodeMap origMap(blockCutTree);
	std::vector<Graph::Node> componentNodes(compCount);
	for (int i = 0; i<compCount; ++i) {
		componentNodes[i] = blockCutTree.addNode();
	}
	for (Graph::NodeIt it(g); it!=lemon::INVALID; ++it) {
		std::vector<int> inComponents;
		for (Graph::IncEdgeIt eIt(g, it); eIt!=lemon::INVALID; ++eIt) {
			int component = compMap[eIt];
			std::vector<Graph::Node>& nodesForComp = origMap[componentNodes[component]];
			if (nodesForComp.empty() || nodesForComp.back()!=it) {
				nodesForComp.push_back(it);
				inComponents.push_back(component);
			}
		}
		if (inComponents.size()>1) {
			Graph::Node cutNode = blockCutTree.addNode();
			origMap[cutNode].push_back(it);
			for (int comp:inComponents) {
				blockCutTree.addEdge(cutNode, componentNodes[comp]);
			}
		}
	}
	return BlockDecomposition(blockCutTree, origMap);
}

CutGenerator::CutStatus SimpleCombCutGen::validate(LinearProgram& lp, const std::vector<double>& solution,
												   CutStatus currentStatus) {
	if (currentStatus!=CutGenerator::valid) {
		return CutStatus::valid;
	}
	Graph fractional;
	//TODO welche brauche ich, welche nicht?
	Graph::NodeMap <city_id> toOrig(fractional);
	std::vector<Graph::Node> toFractional(tsp.getCityCount());
	Graph::EdgeMap <variable_id> edgeToVar(fractional);
	Graph::EdgeMap<double> c(fractional);
	Graph::NodeMap<bool> odd(fractional);
	std::vector<variable_id> oneEdges = tsp_util::createFractionalGraph(tsp, lpData, tolerance, solution,
																		fractional, toOrig,
																		toFractional, edgeToVar, c, odd);
	BlockDecomposition blocks = generateBlocks(fractional);
	std::vector<LinearProgram::Constraint> newConstraints;
	for (Graph::NodeIt mainBlockIt(blocks.blockCutTree); mainBlockIt != lemon::INVALID; ++mainBlockIt) {
		const std::vector<Graph::Node>& mainBlock = blocks[mainBlockIt];
		if (mainBlock.size()>=3) {
			LinearProgram::Constraint constr = checkHandle({mainBlockIt}, fractional, blocks, lp, solution, oneEdges,
														   toOrig, toFractional, odd);
			if (constr.isValid() && constr.isViolated(solution, tolerance)) {
				newConstraints.push_back(constr);
			}
			for (Graph::OutArcIt arcIt1(blocks.blockCutTree, mainBlockIt); arcIt1 != lemon::INVALID; ++arcIt1) {
				Graph::Node commonNode = blocks.blockCutTree.target(arcIt1);
				for (Graph::OutArcIt arcIt2(blocks.blockCutTree, commonNode); arcIt2!=lemon::INVALID; ++arcIt2) {
					Graph::Node otherBlock = blocks.blockCutTree.target(arcIt2);
					if (Graph::id(otherBlock) > Graph::id(mainBlockIt)) {
						const std::vector<Graph::Node>& additionalBlock = blocks[otherBlock];
						if (additionalBlock.size()>=3) {
							constr = checkHandle({mainBlockIt, otherBlock}, fractional, blocks, lp,
												 solution, oneEdges, toOrig, toFractional, odd);
							if (constr.isValid() && constr.isViolated(solution, tolerance)) {
								newConstraints.push_back(constr);
							}
						}
					}
				}
			}
		}
	}
	if (!newConstraints.empty()) {
		lp.addConstraints(newConstraints);
		return CutGenerator::maybe_recalc;
	} else {
		return CutGenerator::valid;
	}
}

LinearProgram::Constraint SimpleCombCutGen::checkHandle
		(const SimpleCombCutGen::Handle& h, const Graph& mainGraph, const SimpleCombCutGen::BlockDecomposition& blocks,
		 LinearProgram& lp, const std::vector<double>& solution, const std::vector<variable_id>& oneEdges,
		 const Graph::NodeMap <city_id>& toTSP, const std::vector<Graph::Node>& toGraph,
		 const Graph::NodeMap<bool>& odd) {
	Graph::NodeMap<bool> inHandle(mainGraph);
	std::vector<Graph::Node> handleNodes;
	for (Graph::Node blockNode:h) {
		for (Graph::Node origNode:blocks[blockNode]) {
			if (!inHandle[origNode]) {
				inHandle[origNode] = true;
				handleNodes.push_back(origNode);
			}
		}
	}
	std::vector<VirtualEdge> teeth = getTeethForHandle(h, oneEdges, blocks, toGraph, toTSP, odd, inHandle, mainGraph,
													   solution);
	if (teeth.size()<3) {
		return LinearProgram::Constraint();
	}
	size_t rhs = handleNodes.size();
	for (const VirtualEdge& e:teeth) {
		rhs += e.size()-1;
	}
	rhs -= (teeth.size()+1)/2;
	double lhs = inducedSum(handleNodes, toTSP, solution);
	for (const VirtualEdge& e:teeth) {
		lhs += inducedSum(e, toTSP, solution);
		if (tolerance.less(rhs, lhs)) {
			break;
		}
	}
	if (!tolerance.less(rhs, lhs)) {
		return LinearProgram::Constraint();
	}
	const size_t none = teeth.size();
	Graph::NodeMap <size_t> inTooth(mainGraph, none);
	for (size_t i = 0; i<teeth.size();) {
		size_t otherRemove = none;
		for (Graph::Node n:teeth[i]) {
			if (inTooth[n]==none && inTooth[n]!=i) {
				inTooth[n] = none;
			} else {
				otherRemove = inTooth[n];
				break;
			}
		}
		if (otherRemove==none) {
			++i;
		} else {
			for (Graph::Node n:teeth[i]) {
				inTooth[n] = none;
				if (!inHandle[n]) {
					inHandle[n] = true;
					handleNodes.push_back(n);
				}
			}
			for (Graph::Node n:teeth[otherRemove]) {
				inTooth[n] = none;
				if (!inHandle[n]) {
					inHandle[n] = true;
					handleNodes.push_back(n);
				}
			}
			teeth.erase(teeth.begin()+otherRemove);
			--i;
			std::swap(teeth[i], teeth.back());
			teeth.pop_back();
		}
	}
	if (teeth.size()>=3) {
		std::vector<city_id> handleTSP;
		handleTSP.reserve(handleNodes.size());
		for (Graph::Node n:handleNodes) {
			handleTSP.push_back(toTSP[n]);
		}
		size_t lpRHS = 0;
		std::vector<variable_id> usedIndices;
		std::vector<double> coeffs(lpData.getVariableCount());
		lpRHS += lpData.sparserInducedSum(handleTSP, coeffs, usedIndices) + 1;
		for (const VirtualEdge& e:teeth) {
			std::vector<city_id> toothTSP;
			toothTSP.reserve(handleNodes.size());
			for (Graph::Node n:handleNodes) {
				toothTSP.push_back(toTSP[n]);
			}
			lpRHS += lpData.sparserInducedSum(toothTSP, coeffs, usedIndices);
		}
		lpRHS -= (1 + teeth.size()) / 2;
		return LinearProgram::Constraint::fromDense(usedIndices, coeffs, LinearProgram::less_eq, lpRHS);
	} else {
		return LinearProgram::Constraint();
	}
}

std::vector<SimpleCombCutGen::VirtualEdge> SimpleCombCutGen::getTeethForHandle
		(const SimpleCombCutGen::Handle& handle, const std::vector<variable_id>& oneEdges,
		 const BlockDecomposition& blocks,
		 const std::vector<Graph::Node>& toGraph, const Graph::NodeMap <city_id>& toCity,
		 const Graph::NodeMap<bool>& odd, const Graph::NodeMap<bool>& inHandle, const Graph& g,
		 const std::vector<double>& solution) {
	std::vector<VirtualEdge> ret;
	//Alle 1-Kanten
	for (variable_id eId:oneEdges) {
		TspLpData::Edge e = lpData.getEdge(eId);
		Graph::Node nodeA = toGraph[e.first];
		Graph::Node nodeB = toGraph[e.second];
		if (inHandle[nodeA]!=inHandle[nodeB]) {
			ret.push_back({nodeA, nodeB});
		}
	}
	VirtualEdge minDiffEdge = {};
	size_t minDiffIndex = 0;
	double minDiffVal = 1;
	Graph::NodeMap<bool> processed(g, false);
	for (Graph::Node component:handle) {
		for (Graph::OutArcIt cutnodeIt(blocks.blockCutTree, component); cutnodeIt!=lemon::INVALID; ++cutnodeIt) {
			assert(blocks[blocks.blockCutTree.target(cutnodeIt)].size()==1);
			Graph::Node u = blocks[blocks.blockCutTree.target(cutnodeIt)][0];
			if (processed[u] || odd[u]) {
				continue;
			}
			processed[u] = true;
			double maxWeight = -1;
			VirtualEdge maxWeightEdge = {};
			city_id uCity = toCity[u];
			//(i) Regular edges (virtual edges with cardinality 2)
			for (Graph::OutArcIt it(g, u); it!=lemon::INVALID; ++it) {
				city_id endCity = toCity[g.target(it)];
				if (!inHandle[g.target(it)]) {
					double cost = solution[lpData.getVariable(uCity, endCity)];
					if (cost>maxWeight) {
						maxWeight = cost;
						maxWeightEdge = {u, g.target(it)};
					}
				}
			}
			for (Graph::OutArcIt blockIt(blocks.blockCutTree, cutnodeIt); blockIt!=lemon::INVALID; ++blockIt) {
				Graph::Node neighborBlock = blocks.blockCutTree.target(blockIt);
				if (std::find(handle.begin(), handle.end(), neighborBlock)==handle.end()) {
					continue;
				}
				const std::vector<Graph::Node>& blockNodes = blocks[neighborBlock];
				//(ii) Blocks of (F, x)
				double cost = inducedSum(blockNodes, toCity, solution)+2-blockNodes.size();
				if (cost>maxWeight) {
					maxWeight = cost;
					maxWeightEdge = blockNodes;
				}
				//(iii) The intersection of a block...
				//TODO das muss auch schöner gehen
				std::vector<Graph::Node> type3 = {u};
				for (Graph::Node pot:blockNodes) {
					for (Graph::OutArcIt neighbors(g, u); neighbors!=lemon::INVALID; ++neighbors) {
						Graph::Node target = g.target(neighbors);
						if (pot==target) {
							type3.push_back(pot);
							break;
						}
					}
				}
				cost = inducedSum(blockNodes, toCity, solution)+2-blockNodes.size();
				if (cost>maxWeight) {
					maxWeight = cost;
					maxWeightEdge = blockNodes;
				}
			}

			if (maxWeight>=0) {
				addAndMinDiff(ret, maxWeightEdge, minDiffEdge, minDiffVal, minDiffIndex, maxWeight);
			}
		}
	}
	if (ret.size()%2==0) {
		if (minDiffVal<0) {
			ret.push_back(minDiffEdge);
		} else {
			if (minDiffIndex>=ret.size()) {
				return {};
			}
			std::swap(ret[minDiffIndex], ret.back());
			ret.pop_back();
		}
	}
	return ret;
}

void SimpleCombCutGen::addAndMinDiff(std::vector<SimpleCombCutGen::VirtualEdge>& out, const VirtualEdge& add,
									 VirtualEdge& minDiffEdge, double& minDiffVal, size_t& minDiffIndex, double weight) {
	double diff = weight-.5;
	if (std::abs(diff)<std::abs(minDiffVal)) {
		minDiffVal = diff;
		minDiffEdge = add;
		minDiffIndex = out.size();
	}
	if (diff>=0) {
		out.push_back(add);
	}
}

//TODO die Graphen sind nicht mehr vollstd, gibt es schönere/schnellere Ansätze?
double SimpleCombCutGen::inducedSum(const std::vector<Graph::Node>& inducing, const Graph::NodeMap <city_id>& toCity,
									const std::vector<double>& solution) {
	double ret = 0;
	for (size_t i1 = 1; i1<inducing.size(); ++i1) {
		city_id city1 = toCity[inducing[i1]];
		for (size_t i2 = 0; i2<i1; ++i2) {
			city_id city2 = toCity[inducing[i2]];
			variable_id var = lpData.getVariable(city1, city2);
			if (var != LinearProgram::invalid_variable) {
				ret += solution[var];
			}
		}
	}
	return ret;
}

SimpleCombCutGen::BlockDecomposition::BlockDecomposition(const Graph& tree,
														 const Graph::NodeMap <std::vector<Graph::Node>>& orig) :
		originalNodes(blockCutTree) {
	lemon::graphCopy(tree, blockCutTree).nodeMap(orig, originalNodes).run();
}

SimpleCombCutGen::BlockDecomposition::BlockDecomposition(const SimpleCombCutGen::BlockDecomposition& old)
		: BlockDecomposition(
		old.blockCutTree, old.originalNodes) {
}
