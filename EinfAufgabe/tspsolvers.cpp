#include <algorithm>
#include <iostream>
#include "tspsolvers.hpp"

using edge_id = Graph::edge_id;
using cost_t = Graph::cost_t;
using node_id = Graph::node_id;
namespace tspsolvers {
	/**
	 * Berechnet eine "relative kurze" Tour in der gegebenen TSP-Instanz, indem immer die kürzeste mögliche Kante
	 * hinzugefügt wird.
	 */
	Tour greedyTSP(const Graph& g) {
		std::vector<edge> sortedEdges(g.getEdgeCount());
		for (edge_id i = 0; i < g.getEdgeCount(); ++i) {
			sortedEdges[i] = g.getEdge(i);
		}
		//Sortieren nach Kosten der entsprechenden Kanten
		std::sort(sortedEdges.begin(), sortedEdges.end(),
				  [](const edge& edgeA, const edge& edgeB) {
					  return edgeA.cost < edgeB.cost;
				  }
		);
		/*
		 * Falls an Knoten i weniger als 2 Kanten anliegen, gibt otherEnd[i] den anderen Knoten in der selben
		 * Zusammenhangskomponente wie i an, der auch Grad <2 hat (Falls i Grad 0 hat, ist dies i selbst), d.h. das
		 * andere Ende des Toursegments. Falls an i 2 Kanten anliegen, ist der Wert beliebig.
		 */
		std::vector<node_id> otherEnd(g.getNodeCount());
		for (node_id i = 0; i < g.getNodeCount(); ++i) {
			otherEnd[i] = i;
		}
		//Speichert die gewählten Kanten, die schon zu einem Knoten inzident sind
		std::vector<std::vector<edge>> edgesAtNode(g.getNodeCount());
		unsigned addedEdges = 0;
		cost_t totalCost = 0;
		for (const edge& e:sortedEdges) {
			//Falls an einem der beiden Enden schon 2 Kanten anliegen, kann die Kante nicht hinzugefügt werden
			std::vector<edge>& edgesAtA = edgesAtNode[e.endA];
			if (edgesAtA.size() >= 2) {
				continue;
			}
			std::vector<edge>& edgesAtB = edgesAtNode[e.endB];
			if (edgesAtB.size() >= 2) {
				continue;
			}
			/*
			 * Wenn die Enden beide Grad <2 haben (also Enden von Toursegmenten sind), können sie genau dann verbunden
			 * werden, wenn sie nicht zum selben Segment gehören.
			 */
			if (otherEnd[e.endA] != e.endB) {
				++addedEdges;
				totalCost += e.cost;
				edgesAtA.push_back(e);
				edgesAtB.push_back(e);
				if (addedEdges == g.getNodeCount() - 1) {
					break;//Tour ist fast vollständig, die letzte Kante ist aber eindeutig bestimmt
				}
				//Enden des neuen Segments setzen
				//newEndA/B zwischenspeichern, falls die Segmente aus einzelnen Knoten bestehen
				node_id newEndA = otherEnd[e.endA];
				node_id newEndB = otherEnd[e.endB];
				otherEnd[newEndA] = newEndB;
				otherEnd[newEndB] = newEndA;
			}
		}
		closeHamiltonPath(g, edgesAtNode);
		std::vector<node_id> order = getOrderFromEdges(edgesAtNode);
		return Tour(order, totalCost, g.getName());
	}

	/**
	 * Bestimmt die Referenz-Tour in der gegebenen TSP-Instanz, bei der die Knoten in der Reihenfolge durchlaufen werden,
	 * in der sie nummeriert sind
	 */
	Tour getReferenceTour(const Graph& g) {
		std::vector<node_id> order;
		order.reserve(g.getNodeCount());
		cost_t totalCost = 0;
		for (node_id curr = 0; curr < g.getNodeCount(); ++curr) {
			order.push_back(curr);
			for (edge_id eid:g.getEdgesAt(curr)) {
				const edge& e = g.getEdge(eid);
				if (e.getOtherNode(curr) == (curr + 1) % g.getNodeCount()) {
					totalCost += e.cost;
					break;
				}
			}
		}
		return Tour(order, totalCost, g.getName());
	}

	/**
	 * Bestimmt die Reihenfolge der Knoten in der Tour, die aus den angegebenen Kanten besteht
	 * @param edgesAtNode: Der Eintrag an der Stelle i speichert die Kanten, die in der Tour zu Knoten i inzident sind
	 */
	std::vector<node_id> getOrderFromEdges(const std::vector<std::vector<edge>>& edgesAtNode) {
		std::vector<node_id> order;
		node_id prevNode = 0;
		node_id currentNode = 0;
		do {
			for (const edge& e:edgesAtNode[currentNode]) {
				node_id otherNode = e.getOtherNode(currentNode);
				if (otherNode != prevNode) {
					order.push_back(currentNode);
					prevNode = currentNode;
					currentNode = otherNode;
					break;
				}
			}
		} while (currentNode != 0);
		return order;
	}

	/**
	 * Fügt die fehlende Kante in einen Hamilton-Pfad in der TSP-Instanz g ein
	 */
	void closeHamiltonPath(const Graph& g, std::vector<std::vector<edge>>& edgesAtNode) {
		node_id firstEnd = g.getNodeCount();
		for (node_id i = 0; i < g.getNodeCount(); ++i) {
			if (edgesAtNode[i].size() == 1) {
				if (firstEnd < g.getNodeCount()) {
					for (edge_id eid:g.getEdgesAt(i)) {
						const edge& e = g.getEdge(eid);
						if (e.getOtherNode(i) == firstEnd) {
							edgesAtNode[i].push_back(e);
							edgesAtNode[firstEnd].push_back(e);
							break;
						}
					}
				} else {
					firstEnd = i;
				}
			}
		}
	}
}