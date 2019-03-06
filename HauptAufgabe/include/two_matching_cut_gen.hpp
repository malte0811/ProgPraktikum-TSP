#ifndef TWO_MATCHING_CUT_GEN_HPP
#define TWO_MATCHING_CUT_GEN_HPP

#include <tsp_instance.hpp>
#include <cut_generator.hpp>
#include <lemon/gomory_hu.h>

class TwoMatchingCutGen : public CutGenerator {
public:
	explicit TwoMatchingCutGen(const TSPInstance& inst, bool contract);

	CutStatus validate(LinearProgram& lp, const std::vector<double>& solution) override;

private:
	using ContractionMap = Graph::NodeMap<std::vector<city_id>>;
	struct XandF {
		std::vector<Graph::Node> x;
		std::vector<Graph::Edge> f;
		double cost;
	};

	void contractPaths(Graph& g, Graph::NodeMap<bool>& odd, Graph::EdgeMap<double>& c,
					   ContractionMap& toOrig);

	void lemma1220(const Graph& graph, std::vector<XandF>& out, const Graph::NodeMap<bool>& odd,
				   const Graph::EdgeMap<double>& c);

	std::vector<Graph::Node> discoverPath(const Graph& graph, Graph::Node start, Graph::Node exclude,
										  const Graph::NodeMap<bool>& odd, Graph::NodeMap<bool>& visited,
										  const Graph::EdgeMap<double>& c, double& firstEdgeVal);

	bool findAndContractPath(Graph& g, Graph::Node start, ContractionMap& toOrig, Graph::NodeMap<bool>& odd,
							 const Graph::EdgeMap<double>& c);

	const TSPInstance& tsp;

	const bool enableContraction;

	const lemon::Tolerance<double> tolerance;
};


#endif
