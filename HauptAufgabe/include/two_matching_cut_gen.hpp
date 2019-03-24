#ifndef TWO_MATCHING_CUT_GEN_HPP
#define TWO_MATCHING_CUT_GEN_HPP

#include <tsp_instance.hpp>
#include <cut_generator.hpp>
#include <lemon/gomory_hu.h>
#include "union_find.hpp"

class TwoMatchingCutGen : public CutGenerator {
public:
	explicit TwoMatchingCutGen(const TSPInstance& inst, bool contract);

	CutStatus validate(LinearProgram& lp, const std::vector<double>& solution) override;

private:
	using ContractionMap = Graph::NodeMap<std::vector<city_id>>;
	struct Blossom {
		std::vector<Graph::Node> x;
		std::vector<Graph::Edge> f;
	};

	void finalizeBlossom(const std::vector<variable_id>& oneEdges, std::vector<bool>& isInX,
						 std::vector<variable_id>& f, size_t& sizeX);

	Blossom calculateAndAddBlossom(const Graph::NodeMap <size_t>& nodeToUF, UnionFind& components,
								   size_t xIndex, const Graph& g, double cutCost,
								   const Graph::EdgeMap<double>& c, const Graph::NodeMap<bool>& odd,
								   const Graph::NodeMap <size_t>& adjacentEDash);

	void contractPaths(Graph& g, Graph::NodeMap<bool>& odd, Graph::EdgeMap<double>& c,
					   ContractionMap& toOrig);

	std::vector<Blossom> lemma1220(const Graph& graph, const Graph::NodeMap<bool>& odd,
								   const Graph::EdgeMap<double>& c);

	std::vector<Graph::Node> discoverPath(const Graph& graph, Graph::Node start,
										  lemon::ListGraphBase::Edge& exclude,
										  const Graph::NodeMap<bool>& odd, Graph::NodeMap<bool>& visited,
										  const Graph::EdgeMap<double>& c);

	bool findAndContractPath(Graph& g, Graph::Node start, ContractionMap& toOrig, Graph::NodeMap<bool>& odd,
							 const lemon::GraphExtender<lemon::ListGraphBase>::EdgeMap<double>& c);

	bool isValidInternalNode(const Graph& g, Graph::Node start, const Graph::EdgeMap<double>& c);

	const TSPInstance& tsp;

	const bool enableContraction;

	const lemon::Tolerance<double> tolerance;
};


#endif
