#ifndef TWO_MATCHING_CUT_GEN_HPP
#define TWO_MATCHING_CUT_GEN_HPP

#include <tsp_instance.hpp>
#include <cut_generator.hpp>
#include <lemon/gomory_hu.h>

class TwoMatchingCutGen : public CutGenerator {
public:
	explicit TwoMatchingCutGen(const TSPInstance& inst);

	bool validate(LinearProgram& lp, const std::vector<double>& solution) override;

private:

	struct XandF {
		explicit XandF(const Graph& graph);
		XandF(const XandF& other);

		const Graph* g;//Needed for the copy-constructor...
		Graph::NodeMap<bool> x;
		Graph::EdgeMap<bool> f;
		double cost;
		city_id sizeX;
		variable_id sizeF;
	};

	//T={}
	std::vector<XandF> lemma1220(const Graph& graph,
								 const Graph::EdgeMap<double>& c,
								 const Graph::EdgeMap<double>& cDash);

	const TSPInstance& tsp;
	Graph workGraph;
	Graph::Snapshot basicState;
	Graph::NodeMap <Graph::Node> origToWork;
	Graph::NodeMap <Graph::Node> workToOrig;
	Graph::Node z;
	Graph::EdgeMap<double> c;
	Graph::EdgeMap<double> cDash;
	lemon::Tolerance<double> tolerance;
};


#endif
