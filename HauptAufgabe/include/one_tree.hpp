#ifndef ONE_TREE_HPP
#define ONE_TREE_HPP


#include <tsp_instance.hpp>
#include <lemon/full_graph.h>
#include <lemon/adaptors.h>
#include "relative_tolerance.hpp"

class OneTree {
public:
	using FullGraph = lemon::FullGraph;
	using Tree = lemon::FilterEdges<FullGraph>;

	explicit OneTree(const TSPInstance& inst);

	void run(CPXENVptr env);

	double getCleanCost() const;

	const std::vector<double>& getPotential() const;

	const Tree& getTree() const;

	const FullGraph& getGraph() const;

	const FullGraph::EdgeMap<bool>& getInTree() const;

	const FullGraph::EdgeMap<double>& getReducedCosts() const;

	double getBound();

private:
	bool generate1Tree();

	void addVariable(LinearProgram& lp);

	void addInitialTrees(LinearProgram& lp);

	const TSPInstance *inst;
	cost_t costOriginal = 0;
	double costPotential = 0;
	double lowerBound;
	//TODO da ist mehr drin als nur das Potential
	std::vector<double> potential;
	FullGraph g;
	FullGraph::EdgeMap<double> costs;
	FullGraph::EdgeMap<bool> inTree;
	FullGraph::EdgeMap<bool> bestInTree;
	std::pair<city_id, city_id> oneNeighbors;
	Tree t;
	RelativeTolerance tolerance;
};


#endif
