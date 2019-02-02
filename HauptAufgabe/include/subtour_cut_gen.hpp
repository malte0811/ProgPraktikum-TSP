#ifndef SUBTOUR_CUT_GEN_HPP
#define SUBTOUR_CUT_GEN_HPP

#include <cut_generator.hpp>
#include <lemon/smart_graph.h>
#include <lemon_fixes/nagamochi_ibaraki.h>
#include <tsp_instance.hpp>

class SubtourCutGen : public CutGenerator {
public:
	explicit SubtourCutGen(const lemon::SmartGraph& graph);

	bool validate(LinearProgram& lp, const std::vector<double>& solution) override;

private:
	const lemon::SmartGraph& origGraph;
	Graph workGraph;
	lemon::SmartGraph::EdgeMap<double> capacity;
	lemon::NagamochiIbaraki<lemon::SmartGraph, lemon::SmartGraph::EdgeMap<double>> minCut;
	lemon::Tolerance<double> tolerance;
};


#endif