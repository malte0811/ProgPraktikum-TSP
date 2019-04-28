#ifndef BLOSSOM_FINDER_HPP
#define BLOSSOM_FINDER_HPP

#include <lemon/tolerance.h>
#include <cstddef>
#include <tsp_instance.hpp>
#include <tsp_utils.hpp>
#include <vector>

using std::size_t;

class UnionFind;

/**
 * Findet verletzte 2-Matching-Constraints/Blüten, falls es solche gibt. Der Algorithmus entspricht grob dem aus Satz
 * 12.21 mit einigen Änderungen:
 * 1. z ist hier nicht nötig, da die Grad-Constraints sind immer erfüllt. Für alle an z anliegenden Kanten e würde
 * c(e)=0 und c'(e)=inf gelten, sie können also nie in der Menge F einer verletzten 2-Matching-Constraint enthalten sein.
 * 2. Der Graph wird vor dem Anwenden von Lemma 12.20 durch Pfad-Kontraktion nach Proposition 4.6 in
 * Padberg, M. & Rinaldi, G. Mathematical Programming (1990) 47: 219. https://doi.org/10.1007/BF01580861
 * 3. Im Algorithmus aus Lemma 12.21 wird hier nicht eine minimale Blüte berechnet, sondern alle aus dem Gomory-Hu-
 * Baum entstehenden Blüten mit Wert echt kleiner 1 (falls die minimalen Blüten Wert kleiner als 1 haben, wird
 * mindestens eine gefunden).
 * 4. Die so gefundenen Blüten werden so verändert, dass F ein Matching ist. Der Algorithmus wird im selben Paper wie
 * oben ohne Beweis gegeben. TODO: Quelle mit Beweis finden oder selbst einen schreiben
 */
class BlossomFinder {
public:
	struct Blossom {
		std::vector<Graph::Node> handle;
		std::vector<Graph::Edge> teeth;

		void replaceByMapped(const tsp_util::ContractionMapGraph& nodeMap, const Graph::EdgeMap <Graph::Edge>& edgeMap);

		bool isProperBlossom() const;
	};

	BlossomFinder(const Graph& g, Graph::EdgeMap<double>& capacities, lemon::Tolerance<double> tolerance,
				  bool shouldContractPaths);

	std::vector<Blossom> findViolatedBlossoms();

private:

	void setupFractionalGraph();

	bool finalizeBlossom(Blossom& b);

	Blossom calculateBlossomFor(const Graph::NodeMap <size_t>& nodeToUF, UnionFind& components, size_t xIndex,
								double cutCost, const Graph::NodeMap <size_t>& adjacentEDash);

	void contractPaths();

	std::vector<Blossom> lemma1220();

	std::vector<Graph::Node> discoverPath(Graph::Node start, Graph::Edge& exclude, Graph::NodeMap<bool>& visited);

	bool findAndContractPath(Graph::Node start);

	bool isValidInternalNode(Graph::Node start);

	lemon::Tolerance<double> tolerance;
	//Der Eingabegraph, in dem Blüten gefunden werden sollen
	const Graph& inputGraph;
	/*
	 * Der fraktionale Graph (d.h. der Eingabegraph ohne Kanten mit Wert 1). Kanten mit Wert 0 kommen in den Eingabe-
	 * graphen nicht vor (würden aber kein Problem darstellen). Wird während des Algorithmus kontrahiert, falls
	 * shouldContractPaths wahr ist
	 */
	Graph fractionalGraph;
	//Kantenkapazitäten auf dem Eingabegraphen
	Graph::EdgeMap<double>& capacitiesMain;
	//Kantenkapazitäten auf dem fraktionalen Graphen
	Graph::EdgeMap<double> capacitiesFractional;
	//Ungerade Knoten im fraktionalen Graph (d.h. solche, an denen eine ungerade Anzahl an Kanten mit Wert 1 entfernt wurden)
	Graph::NodeMap<bool> oddNodes;
	/*
	 * Ordnet den (möglicherweise kontrahierten) Knoten im fraktionalen Graphen die entsprechenden Knoten im
	 * Eingabegraphen zu
	 */
	tsp_util::ContractionMapGraph contraction;
	//Ordnet den Kanten im fraktionalen Graphen die entsprechenden Kanten im Eingabegraphen zu
	Graph::EdgeMap <Graph::Edge> fractToMainEdges;
	//Die Kanten mit Wert 1, die im Eingabegraphen, aber nicht im fraktionalen Graphen enthalten sind
	std::vector<Graph::Edge> oneEdges;
	//Anzahl der Knoten im Eingabe- bzw. im fraktionalen Graph
	const size_t nodeCount;
};


#endif