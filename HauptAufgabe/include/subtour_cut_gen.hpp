#ifndef SUBTOUR_CUT_GEN_HPP
#define SUBTOUR_CUT_GEN_HPP

#include <cut_generator.hpp>
#include <lemon/smart_graph.h>
#include <lemon_fixes/nagamochi_ibaraki.h>
#include <tsp_instance.hpp>

class SubtourCutGen : public CutGenerator {
public:
	explicit SubtourCutGen(const TSPInstance& inst);
	bool validate(LinearProgram& lp, const std::vector<double>& solution) override;

private:
	void addConnectivityConstraints(LinearProgram& lp);

	void addCutConstraint(LinearProgram& lp);
	const TSPInstance& tsp;
	Graph workGraph;
	Graph::Snapshot baseState;
	Graph::NodeMap <Graph::Node> origToWork;
	Graph::NodeMap <Graph::Node> workToOrig;
	Graph::EdgeMap<double> capacity;
	lemon::NagamochiIbaraki<Graph, Graph::EdgeMap<double>> minCut;
	lemon::Tolerance<double> tolerance;
};


#endif