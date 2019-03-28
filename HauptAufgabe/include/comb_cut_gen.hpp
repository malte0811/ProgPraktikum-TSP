#ifndef COMB_CUT_GEN_HPP
#define COMB_CUT_GEN_HPP

#include <cut_generator.hpp>
#include <tsp_instance.hpp>
#include <lemon/connectivity.h>

class CombCutGen : public CutGenerator {
public:
	explicit CombCutGen(const TSPInstance& inst);

	CutStatus validate(LinearProgram& lp, const std::vector<double>& solution, CutStatus currentStatus) override;

private:
	struct BlockDecomposition {
		using OriginalNodeMap = Graph::NodeMap<std::vector<Graph::Node>>;

		BlockDecomposition(const BlockDecomposition& old);

		BlockDecomposition(const Graph& tree, const Graph::NodeMap <std::vector<Graph::Node>>& orig);

		const std::vector<Graph::Node>& operator[](const Graph::Node& n) const {
			return originalNodes[n];
		}

		Graph blockCutTree;
		//1 Element: Cutnode, >1 Element: Komponente
		OriginalNodeMap originalNodes;
	};
	//Knoten im blockCutTree
	using Handle = std::vector<Graph::Node>;
	//Knoten im fraktionalen Graphen
	using VirtualEdge = std::vector<Graph::Node>;

	BlockDecomposition generateBlocks(const Graph& g);

	LinearProgram::Constraint checkHandle(const CombCutGen::Handle& h, const Graph& mainGraph,
										  const CombCutGen::BlockDecomposition& blocks, LinearProgram& lp,
										  const std::vector<double>& solution,
										  const std::vector<variable_id>& oneEdges,
										  const Graph::NodeMap <city_id>& toTSP,
										  const std::vector<Graph::Node>& toGraph, const Graph::NodeMap<bool>& odd);

	const TSPInstance& tsp;
	const lemon::Tolerance<double> tolerance;

	std::vector<CombCutGen::VirtualEdge> getTeethForHandle(const CombCutGen::Handle& handle,
														   const std::vector<variable_id>& oneEdges,
														   const BlockDecomposition& blocks,
														   const std::vector<Graph::Node>& toGraph,
														   const Graph::NodeMap <city_id>& toCity,
														   const Graph::NodeMap<bool>& odd,
														   const Graph::NodeMap<bool>& inHandle, const Graph& g,
														   const std::vector<double>& solution);

	double inducedSum(const std::vector<Graph::Node>& inducing, const Graph::NodeMap <city_id>& toCity,
					  const std::vector<double>& solution);

	void inducedSum(const std::vector<Graph::Node>& inducing, const Graph::NodeMap <city_id>& toCity,
					std::vector<variable_id>& out);

	void addAndMinDiff(std::vector<CombCutGen::VirtualEdge>& out, const VirtualEdge& add,
					   VirtualEdge& minDiffEdge, double& minDiffVal, size_t& minDiffIndex, double weight);
};


#endif