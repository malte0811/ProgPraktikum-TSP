#ifndef COMB_CUT_GEN_HPP
#define COMB_CUT_GEN_HPP

#include <contraction_rule.hpp>

class CombCutGen {
public:
	CombCutGen() :
			onePath(contractOnePath), triangle2(contractTriangle2), square3(contractSquare3),
			oneSquare(contractOneSquare), triangleGE05(contractTriangleGE05) {}

private:
	ContractionRule onePath;
	ContractionRule triangle2;
	ContractionRule square3;
	ContractionRule oneSquare;
	ContractionRule triangleGE05;

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