#ifndef TWO_MATCHING_CUT_GEN_HPP
#define TWO_MATCHING_CUT_GEN_HPP

#include <tsp_instance.hpp>
#include <cut_generator.hpp>
#include <lemon/gomory_hu.h>

class TwoMatchingCutGen : public CutGenerator {
public:
	explicit TwoMatchingCutGen(const TSPInstance& inst, bool contract);

	bool validate(LinearProgram& lp, const std::vector<double>& solution) override;

private:

	struct XandF {
		std::vector<Graph::Node> x;
		std::vector<Graph::Edge> f;
		double cost;
	};

	void contractPaths(Graph& g, Graph::NodeMap<bool>& odd, Graph::EdgeMap<double>& c,
					   Graph::EdgeMap <variable_id>& vars,
					   Graph::NodeMap <std::vector<Graph::Node>>& toOrig);

	void lemma1220(const Graph& graph, std::vector<XandF>& out, const Graph::NodeMap<bool>& odd,
				   const Graph::EdgeMap<double>& c);

	void discoverPath(const Graph& graph, Graph::Node start, Graph::Node exclude,
					  const Graph::NodeMap<bool>& odd, Graph::NodeMap<bool>& visited,
					  std::vector<Graph::Node>& path, const Graph::EdgeMap<double>& c,
					  const Graph::NodeMap <std::vector<Graph::Node>>& toOrig, double& firstEdgeVal);
	const TSPInstance& tsp;

	const lemon::Tolerance<double> tolerance;

	const bool enableContraction;
};


#endif
