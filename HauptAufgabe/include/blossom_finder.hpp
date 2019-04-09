#ifndef BLOSSOM_FINDER_HPP
#define BLOSSOM_FINDER_HPP


#include <vector>
#include <tsp_instance.hpp>
#include <union_find.hpp>
#include <tsp_utils.hpp>

class BlossomFinder {
public:
	struct Blossom {
		std::vector<Graph::Node> handle;
		std::vector<Graph::Edge> teeth;

		void replaceByMapped(const tsp_util::ContractionMapGraph& nodeMap, const Graph::EdgeMap <Graph::Edge>& edgeMap);

		bool isProperBlossom() const;
	};

	BlossomFinder(const Graph& g, Graph::EdgeMap<double>& capacities, lemon::Tolerance<double> tolerance,
				  bool contractPaths);

	std::vector<Blossom> findViolatedBlossoms();

private:

	void finalizeBlossom(Blossom& b, const std::vector<Graph::Edge>& oneEdges);

	Blossom calculateAndAddBlossom(const Graph::NodeMap <size_t>& nodeToUF, UnionFind& components, size_t xIndex,
								   double cutCost, const Graph::NodeMap<bool>& odd,
								   const Graph::NodeMap <size_t>& adjacentEDash);

	void contractPaths(Graph::NodeMap<bool>& odd, tsp_util::ContractionMapGraph& toOrig);

	std::vector<Blossom> lemma1220(const Graph::NodeMap<bool>& odd);

	std::vector<Graph::Node> discoverPath(Graph::Node start, Graph::Edge& exclude, const Graph::NodeMap<bool>& odd,
										  Graph::NodeMap<bool>& visited);

	bool findAndContractPath(Graph::Node start, tsp_util::ContractionMapGraph& toOrig, Graph::NodeMap<bool>& odd);

	bool isValidInternalNode(const Graph& g, Graph::Node start, const Graph::EdgeMap<double>& c);

	lemon::Tolerance<double> tolerance;
	const Graph& mainGraph;
	Graph fractionalGraph;
	Graph::EdgeMap<double>& capacitiesMain;
	Graph::EdgeMap<double> capacitiesFractional;
	const bool shouldContractPaths;
	const size_t nodeCount;
};


#endif