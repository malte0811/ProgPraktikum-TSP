#ifndef COMB_CUT_GEN_HPP
#define COMB_CUT_GEN_HPP

#include <contraction_rule.hpp>
#include <comb_heuristic.hpp>

/*
 * Findet heuristisch verletzte Kamm-Ungleichungen mit der Heuristik aus: Martin Grötschel, Olaf Holland: "Solution of
 * large-scale symmetric travelling salesman problems", 1991
 */
class CombCutGen : public CutGenerator {
public:
	CombCutGen(const TSPInstance& tsp, const TspLpData& lpData);

	CutStatus validate(LinearProgram& lp, const std::vector<double>& solution, CutStatus currentStatus) override;

private:
	//Die Transformationen und die daraus entstehenden Heuristiken aus dem Paper
	ContractionRule onePath;
	ContractionRule triangle2;
	ContractionRule square3;
	ContractionRule oneSquare;
	ContractionRule triangleGE05;
	CombHeuristic heuristic1;
	CombHeuristic heuristic2;
	CombHeuristic heuristic3;
	CombHeuristic heuristic4;
	const TSPInstance& tsp;
	const TspLpData& lpData;
	const lemon::Tolerance<double> tolerance;

	LinearProgram::Constraint getContraintFor(const CombHeuristic::Comb& c);

	static ContractionRule::Contraction findOnePath(const Graph& g, const std::vector<Graph::Node>& possibleNodes,
													const std::vector<Graph::Edge>&,
													const Graph::EdgeMap<double>& costs,
													const Graph::NodeMap<bool>& used);

	static ContractionRule::Contraction findTriangle2(const Graph& g, const std::vector<Graph::Node>&,
													  const std::vector<Graph::Edge>& possibleOneEdges,
													  const Graph::EdgeMap<double>& costs,
													  const Graph::NodeMap<bool>&);

	static ContractionRule::Contraction findSquare3(const Graph& g, const std::vector<Graph::Node>& possibleNodes,
													const std::vector<Graph::Edge>&,
													const Graph::EdgeMap<double>& costs,
													const Graph::NodeMap<bool>& used);

	static ContractionRule::Contraction findOneSquare(const Graph& g, const std::vector<Graph::Node>&,
													  const std::vector<Graph::Edge>& possibleOneEdges,
													  const Graph::EdgeMap<double>& costs, const Graph::NodeMap<bool>&);

	static ContractionRule::Contraction findTriangleGE05(const Graph& g, const std::vector<Graph::Node>&,
														 const std::vector<Graph::Edge>& possibleOneEdges,
														 const Graph::EdgeMap<double>& costs,
														 const Graph::NodeMap<bool>&);
};


#endif