#ifndef COMB_CUT_GEN_HPP
#define COMB_CUT_GEN_HPP

#include <contraction_rule.hpp>
#include <comb_heuristic.hpp>

class CombCutGen : public CutGenerator {
public:
	CombCutGen(const TSPInstance& tsp, const TspLpData& lpData);

	CutStatus validate(LinearProgram& lp, const std::vector<double>& solution, CutStatus currentStatus) override;

private:
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

	LinearProgram::Constraint getContraintFor(const CombHeuristic::Comb& c);

	static ContractionRule::Contraction contractOnePath(const Graph& g, const std::vector<Graph::Node>& possibleNodes,
														const std::vector<Graph::Edge>& possibleOneEdges,
														const Graph::EdgeMap<double>& costs,
														const Graph::NodeMap<bool>& used);

	static ContractionRule::Contraction contractTriangle2(const Graph& g,
														  const std::vector<Graph::Node>& possibleNodes,
														  const std::vector<Graph::Edge>& possibleOneEdges,
														  const Graph::EdgeMap<double>& costs,
														  const Graph::NodeMap<bool>& used);

	static ContractionRule::Contraction contractSquare3(const Graph& g, const std::vector<Graph::Node>& possibleNodes,
														const std::vector<Graph::Edge>& possibleOneEdges,
														const Graph::EdgeMap<double>& costs,
														const Graph::NodeMap<bool>& used);

	static ContractionRule::Contraction contractOneSquare(const Graph& g, const std::vector<Graph::Node>& possibleNodes,
														  const std::vector<Graph::Edge>& possibleOneEdges,
														  const Graph::EdgeMap<double>& costs,
														  const Graph::NodeMap<bool>& used);

	static ContractionRule::Contraction contractTriangleGE05(const Graph& g,
															 const std::vector<Graph::Node>& possibleNodes,
															 const std::vector<Graph::Edge>& possibleOneEdges,
															 const Graph::EdgeMap<double>& costs,
															 const Graph::NodeMap<bool>& used);
};


#endif