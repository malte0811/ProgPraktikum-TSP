#ifndef COMB_CUT_GEN_HPP
#define COMB_CUT_GEN_HPP

#include <cut_generator.hpp>
#include <tsp_instance.hpp>
#include <lemon/connectivity.h>
#include "tsp_lp_data.hpp"

class SimpleCombCutGen : public CutGenerator {
public:
	explicit SimpleCombCutGen(const TSPInstance& inst, const TspLpData& lpData);

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

	const TSPInstance& tsp;
	const TspLpData& lpData;
	const lemon::Tolerance<double> tolerance;

	BlockDecomposition generateBlocks(const Graph& g);

	LinearProgram::Constraint checkHandle(const SimpleCombCutGen::Handle& h, const Graph& mainGraph,
										  const SimpleCombCutGen::BlockDecomposition& blocks, LinearProgram& lp,
										  const std::vector<double>& solution,
										  const std::vector<variable_id>& oneEdges,
										  const Graph::NodeMap <city_id>& toTSP,
										  const std::vector<Graph::Node>& toGraph, const Graph::NodeMap<bool>& odd);

	std::vector<SimpleCombCutGen::VirtualEdge> getTeethForHandle(const SimpleCombCutGen::Handle& handle,
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

	void addAndMinDiff(std::vector<SimpleCombCutGen::VirtualEdge>& out, const VirtualEdge& add,
					   VirtualEdge& minDiffEdge, double& minDiffVal, size_t& minDiffIndex, double weight);
};


#endif