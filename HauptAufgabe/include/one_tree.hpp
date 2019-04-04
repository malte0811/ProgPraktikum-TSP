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

	void run();

	double getCleanCost() const;

	const std::vector<double>& getPotential() const;

	const Tree& getTree() const;

	const FullGraph& getGraph() const;

	const FullGraph::EdgeMap<bool>& getInTree() const;

	const FullGraph::EdgeMap<double>& getReducedCosts() const;

	double getBound();

	std::pair<city_id, city_id> getOneNeighbors() const;

private:
	using Sequence = std::vector<std::pair<FullGraph::Edge, double>>;

	void generate1Tree();

	const TSPInstance *inst;
	cost_t costOriginal = 0;
	double costPotential = 0;
	double lowerBound;
	Sequence sortedEdges;

	std::vector<double> potential;
	std::vector<int> lastV;
	FullGraph g;
	FullGraph::EdgeMap<double> costs;
	FullGraph::EdgeMap<bool> inTree;
	FullGraph::EdgeMap<bool> bestInTree;
	std::pair<city_id, city_id> oneNeighbors;
	Tree t;
	RelativeTolerance tolerance;

	bool updatePotential(double stepSize);
};


#endif
