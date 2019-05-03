#include <blossom_finder.hpp>
#include <cassert>
#include <lemon/core.h>
#include <lemon/gomory_hu.h>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <utility>
#include <tsp_instance.hpp>
#include <tsp_utils.hpp>
#include <union_find.hpp>

BlossomFinder::BlossomFinder(const Graph& g, Graph::EdgeMap<double>& capacities, lemon::Tolerance<double> tolerance,
							 bool shouldContractPaths) :
		tolerance(tolerance), inputGraph(g), capacitiesMain(capacities), capacitiesFractional(fractionalGraph),
		oddNodes(fractionalGraph, false), contraction(fractionalGraph),
		fractToMainEdges(fractionalGraph, lemon::INVALID), nodeCount(lemon::countNodes(inputGraph)) {
	setupFractionalGraph();
	if (shouldContractPaths) {
		contractPaths();
	}
}

/**
 * Berechnet möglichst viele von den gegebenen Gewichten verletzte Blüten. Sollte nur einmal aufgerufen werden.
 * @return Ein Vector von verletzten Blüten
 */
std::vector<BlossomFinder::Blossom> BlossomFinder::findViolatedBlossoms() {
	std::vector<Blossom> violated = lemma1220();
	for (size_t i = 0; i < violated.size();) {
		Blossom& b = violated[i];
		b.replaceByMapped(contraction, fractToMainEdges);
		if (finalizeBlossom(b)) {
			++i;
			assert(b.teeth.size() % 2 == 1);
		} else {
			b = std::move(violated.back());
			violated.pop_back();
		}
	}
	return violated;
}

/**
 * Berechnet wie in Lemma 12.20 Blüten mit Wert kleiner als 1. c'(e) ist für alle Kanten 1-c(e)
 * @return ein Vector mit Blüten mit Wert kleiner als 1. Falls es solche gibt, ist der Vector nicht leer.
 */
std::vector<BlossomFinder::Blossom> BlossomFinder::lemma1220() {
	//d wie im Beweis von Lemma 12.20
	Graph::EdgeMap<double> d(fractionalGraph);
	//Die Anzahl der inzidenten Kanten in E'
	Graph::NodeMap <size_t> adjacentEDash(fractionalGraph);
	//d und adjacentEDash berechnen
	for (Graph::EdgeIt it(fractionalGraph); it != lemon::INVALID; ++it) {
		if (capacitiesFractional[it] > 0.5) {
			++adjacentEDash[fractionalGraph.v(it)];
			++adjacentEDash[fractionalGraph.u(it)];
			d[it] = 1 - capacitiesFractional[it];
		} else {
			d[it] = capacitiesFractional[it];
		}
		assert(d[it] > -0.01);
	}
	//Gomory-Hu-Baum/Arboreszenz bezüglich d berechnen
	lemon::GomoryHu<Graph, Graph::EdgeMap < double>>
	gh(fractionalGraph, d);
	gh.run();
	/*
	 * Die Arboreszenz durchlaufen. Eine einfachere Implementierung, bei der wiederholt minCutMap aufgerufen wird, ist
	 * deutlich langsamer. Die Arboreszenz wird von den Blättern aus zur Wurzel durchlaufen. Die Knoten auf der "unteren"
	 * Seite der aktuellen Kante sind als Menge in components gespeichert, nach dem Bearbeiten einer Kante werden die
	 * Komponenten an den Enden der Kante zusammengefügt.
	 */
	//Die Anzahl der Kinder eines Knotens in der Gomory-Hu-Arboreszenz
	Graph::NodeMap <size_t> childCount(fractionalGraph);
	//Die Anzahl der Knoten im Graphen bzw. in der Arboreszenz
	for (Graph::NodeIt it(fractionalGraph); it != lemon::INVALID; ++it) {
		Graph::Node pred = gh.predNode(it);
		if (pred != lemon::INVALID) {
			++childCount[pred];
		}
	}
	//Ordnet jedem Knoten eine ID im Union-Find zu
	Graph::NodeMap <size_t> nodeToUF(fractionalGraph);
	UnionFind components(nodeCount);
	/*
	 * Die Blätter der Arboreszenz. Wenn eine Kante behandelt wurde, wird ihr "unterer" Endknoten aus der Arboreszenz
	 * "entfernt", der "obere" Endknoten wird also evtl zu einem Blatt und wird zu diesem Vector hinzugefügt.
	 */
	std::vector<Graph::Node> leaves;
	{
		size_t nextIndex = 0;
		for (Graph::NodeIt it(fractionalGraph); it != lemon::INVALID; ++it, ++nextIndex) {
			if (childCount[it] == 0) {
				leaves.push_back(it);
			}
			nodeToUF[it] = nextIndex;
		}
	}
	std::vector<Blossom> ret;
	while (!leaves.empty()) {
		Graph::Node currentNode = leaves.back();
		leaves.pop_back();
		Graph::Node pred = gh.predNode(currentNode);
		//Für die Wurzel der Arboreszenz gibt es keinen Vorgänger und auch keine zu betrachtende Kante
		if (pred == lemon::INVALID) {
			continue;
		}
		//Repräsentant der "unteren" Seite der aktuellen Kante
		const size_t xIndex = components.find(nodeToUF[currentNode]);
		double cutCost = gh.predValue(currentNode);
		//Der Wert der Blüte ist mindestens der Wert des Schnitts
		if (tolerance.less(cutCost, 1)) {
			Blossom b = calculateBlossomFor(nodeToUF, components, xIndex, cutCost, adjacentEDash);
			if (!b.handle.empty()) {
				ret.push_back(b);
			}
		}
		components.mergeRoots(xIndex, components.find(nodeToUF[pred]));
		--childCount[pred];
		if (childCount[pred] == 0) {
			leaves.push_back(pred);
		}
	}
	return ret;
}

/**
 * Kontrahiert so lange alternierende Pfade in G, bis es keine mehr gibt, und speichert die Kontraktionen in toOrig
 */
void BlossomFinder::contractPaths() {
	bool contracted;
	do {
		contracted = false;
		for (Graph::NodeIt start(fractionalGraph); start != lemon::INVALID && !contracted; ++start) {
			if (oddNodes[start]) {
				contracted = findAndContractPath(start);
			}
		}
	} while (contracted);
}

/**
 * Findet und kontrahiert einen start enthaltenden alternierenden Pfad in G, falls es einen solchen gibt
 * @return true, genau dann wenn ein Pfad kontrahiert wurde
 */
bool BlossomFinder::findAndContractPath(Graph::Node start) {
	if (!isValidInternalNode(start)) {
		return false;
	}
	Graph::NodeMap<bool> inPath(fractionalGraph);
	//Wird im ersten discoverPath-Aufruf auf die erste Kante des linken Teilpfades gesetzt
	Graph::Edge firstEdge = lemon::INVALID;
	//Teilpfad in eine Richtung ("links") finden
	std::vector<Graph::Node> left = discoverPath(start, firstEdge, inPath);
	std::vector<Graph::Node> right;
	//Ein isolierter Kreis aus ungeraden Knoten
	bool isCycle = left.front() == left.back();
	//Falls der Pfad noch kein Kreis ist, gibt es einen rechten Teilpfad
	if (!isCycle) {
		//Anderen Teilpfad finden
		right = discoverPath(start, firstEdge, inPath);
		//Kreis mit einem geraden Knoten
		isCycle = right.back() == left.back();
	}
	//Der Knoten, der von der kontrahierten Menge übrig bleibt
	Graph::Node remainingNode = left.back();
	//Die restlichen Knoten der kontrahierten Menge:
	//Der letzte Knoten von left wird nie kontrahiert (er bleibt übrig), der erste Knoten ist entweder auch der letzte
	//Knoten oder gleich dem ersten Knoten der rechten Seite
	std::vector<Graph::Node> toRemove(left.begin() + 1, left.end() - 1);
	if (right.size() > 1) {
		//Der letzte rechte Knoten wird nicht hinzugefügt: Falls der Pfad ein Kreis ist, ist dies der verbleibende Knoten;
		//sonst kann nur P-u_1 kontrahiert werden
		toRemove.insert(toRemove.end(), right.begin(), right.end() - 1);
	}
	if (!isCycle) {
		//Falls der Pfad kein Kreis ist, muss die Kante zwischen dem letzten kontrahierten Knoten und u_1 verbleiben
		//Es kann nicht einfach die Kante gelöscht und durch eine neue ersetzt werden, da sonst Werte in EdgeMap's
		//verloren gehen (z.B. toVariable in validate)
		Graph::Node lastContracted = right[right.size() - 2];
		Graph::Node otherRemaining = right.back();
		for (Graph::OutArcIt it(fractionalGraph, lastContracted); it != lemon::INVALID; ++it) {
			Graph::Node target = fractionalGraph.target(it);
			if (target == otherRemaining) {
				if (fractionalGraph.u(it) == lastContracted) {
					fractionalGraph.changeU(it, remainingNode);
				} else {
					fractionalGraph.changeV(it, remainingNode);
				}
				break;
			}
		}
	}
	//Ob der durch die Kontraktion entstehende Knoten ungerade ist
	bool resultOdd = oddNodes[remainingNode];
	std::vector<Graph::Node>& contractedSet = contraction[remainingNode];
	//Kontrahieren
	for (Graph::Node remove:toRemove) {
		assert(remove != remainingNode);
		if (oddNodes[remove]) {
			resultOdd = !resultOdd;
		}
		contractedSet.insert(contractedSet.end(), contraction[remove].begin(), contraction[remove].end());
		fractionalGraph.erase(remove);
	}
	oddNodes[remainingNode] = resultOdd;
	return true;
}

/**
 * Stellt fest, ob ein Knoten als innerer Knoten eines alternierenden Pfades genutzt werden kann
 * @param g Der Graph, in dem der alternierende Pfad existieren soll
 * @param start der Knoten, der im Pfad enthalten sein soll
 * @param c Die Kantengewichte
 * @return wahr genau dann, wenn der Knoten als innerer Knoten in Frage kommt
 */
bool BlossomFinder::isValidInternalNode(Graph::Node start) {
	if (!oddNodes[start]) {
		return false;
	}
	size_t degree = 0;
	double edgeSum = 0;
	for (Graph::OutArcIt it(fractionalGraph, start);
		 it != lemon::INVALID && degree <= 2 && !tolerance.less(1, edgeSum); ++it) {
		++degree;
		edgeSum += capacitiesFractional[it];
	}
	return degree == 2 && !tolerance.different(edgeSum, 1);
}

/**
 * Findet einen Teil eines alternierenden Pfades
 * @param start der Ausgangsknoten. Muss als innerer Knoten eines alternierenden Pfades gültig sein.
 * @param exclude Entweder lemon::INVALID: in diesem Fall wird die erste Kante im Teilpfad hier gespeichert; oder eine
 * Kante, die nicht im Teilpfad genutzt werden darf
 * @param visited Gibt an, ob ein Knoten schon im Pfad enthalten ist
 * @return den Teilpfad eines alternierenden Pfades mit den angegebenen Eigenschaften. Das erste Element ist immer start
 */
std::vector<Graph::Node> BlossomFinder::discoverPath(Graph::Node start, Graph::Edge& exclude,
													 Graph::NodeMap<bool>& visited) {
	//Die Knoten auf dem Pfad
	std::vector<Graph::Node> path;
	Graph::Node current = start;
	//Die letzte Kante im Pfad, also die Kante, die nicht als nächste Kante genutzt werden kann.
	Graph::Edge lastEdge = exclude;
	//Gibt an, ob der aktuelle Knoten ein innerer Knoten in einem alternierenden Pfad sein kann.
	bool canContinuePath;
	do {
		path.push_back(current);
		visited[current] = true;
		Graph::Node nextNode = lemon::INVALID;
		for (Graph::OutArcIt it(fractionalGraph, current); it != lemon::INVALID && nextNode == lemon::INVALID; ++it) {
			//current ist ein gültiger innerer Knoten, hat also Grad 2 und Kantensumme 1. it ist nicht die vorherige
			//Kante im Pfad, also muss it die nächste Kante des Pfades sein
			if (lastEdge != it) {
				nextNode = fractionalGraph.target(it);
				lastEdge = it;
				//Erster Fall der Definition von exclude
				if (exclude == lemon::INVALID) {
					exclude = it;
				}
			}
		}
		//Prüfen, ob der nächste Knoten ein gültiger innerer Knoten ist
		assert(fractionalGraph.valid(nextNode));
		canContinuePath = fractionalGraph.valid(nextNode) && !visited[nextNode] && isValidInternalNode(nextNode);
		current = nextNode;
	} while (canContinuePath);
	if (current != lemon::INVALID) {
		path.push_back(current);
		visited[current] = true;
	}
	return path;
}

/**
 * Fügt die 1-Kanten zur Blüte hinzu und wandelt sie in eine Blüte mit selbem Wert um, deren Zähne ein Matching bilden.
 * Der Rückgabewert ist wahr, falls eine gültige Blüte entstanden ist.
 */
bool BlossomFinder::finalizeBlossom(Blossom& b) {
	const size_t minHandleSize = 3;
	std::vector<Graph::Edge>& teeth = b.teeth;
	std::vector<Graph::Node>& handle = b.handle;
	//Falls der Griff zu klein oder zu groß ist, ist die Blüte ungültig
	if (handle.size() < minHandleSize || handle.size() > nodeCount - minHandleSize) {
		return false;
	}
	Graph::NodeMap<int> handleIndex(inputGraph, -1);
	for (size_t i = 0; i < handle.size(); ++i) {
		handleIndex[handle[i]] = i;
	}
	//Alle 1-Kanten im Schnitt vom Griff sind Zähne
	for (Graph::Edge e:oneEdges) {
		bool uInHandle = handleIndex[inputGraph.u(e)] >= 0;
		bool vInHandle = handleIndex[inputGraph.v(e)] >= 0;
		if (uInHandle != vInHandle) {
			teeth.push_back(e);
		}
	}
	assert(teeth.size() % 2 == 1);
	//Gibt an, ob eine Variable/Kante im TSP-Graphen ein Zahn ist
	Graph::EdgeMap<bool> isTooth(inputGraph, false);
	{
		//Anzahl der zu einem Knoten inzidenten Zähne
		Graph::NodeMap <size_t> incident(inputGraph);
		/*
		 * Setzt die Einträge von isTooth und stellt sicher, dass zu jedem Knoten maximal 2 Zähne inzident sind. Wenn
		 * mehr als 3 Zähne inzident sind, kann die Blüte nicht verletzt sein, durch Rundungsfehler werden trotzdem
		 * manchmal solche Blüten berechnet.
		 */
		for (Graph::Edge tooth:teeth) {
			isTooth[tooth] = true;
			++incident[inputGraph.u(tooth)];
			if (incident[inputGraph.u(tooth)] >= 3) {
				return false;
			}
			++incident[inputGraph.v(tooth)];
			if (incident[inputGraph.v(tooth)] >= 3) {
				return false;
			}
		}
	}
	//Ordnet jedem Knoten den inzidenten Zahn zu
	Graph::NodeMap <Graph::Edge> incidentTooth(inputGraph, lemon::INVALID);
	for (Graph::Edge e:teeth) {
		for (Graph::Node end:{inputGraph.u(e), inputGraph.v(e)}) {
			if (incidentTooth[end] == lemon::INVALID) {
				//Es gibt noch keinen Zahn, der zu end inzident ist
				incidentTooth[end] = e;
			} else {
				//Die Zähne entfernen, die ein gemeinsames Ende haben
				Graph::Edge oldEdge = incidentTooth[end];
				incidentTooth[inputGraph.u(oldEdge)] = lemon::INVALID;
				incidentTooth[inputGraph.v(oldEdge)] = lemon::INVALID;
				incidentTooth[inputGraph.u(e)] = lemon::INVALID;
				incidentTooth[inputGraph.v(e)] = lemon::INVALID;
				isTooth[oldEdge] = false;
				isTooth[e] = false;
				//Das gemeinsame Ende aus dem Griff entfernen bzw zum Griff hinzufügen
				if (handleIndex[end] >= 0) {
					assert(handle[handleIndex[end]] == end);
					handleIndex[handle.back()] = handleIndex[end];
					std::swap(handle[handleIndex[end]], handle.back());
					handle.pop_back();
					handleIndex[end] = -1;
				} else {
					handleIndex[end] = handle.size();
					handle.push_back(end);
				}
				break;
			}
		}
	}
	//Zähne, die zu Erzeugen des Matchings entfernt wurden, auch aus den in der Blüte gespeicherten Zähnen entfernen
	size_t i = 0;
	while (i < teeth.size()) {
		Graph::Edge tooth = teeth[i];
		if (!isTooth[tooth]) {
			std::swap(teeth[i], teeth.back());
			teeth.pop_back();
		} else {
			assert((handleIndex[inputGraph.u(tooth)] >= 0) != (handleIndex[inputGraph.v(tooth)] >= 0));
			++i;
		}
	}
	return handle.size() >= minHandleSize && handle.size() <= nodeCount - minHandleSize;
}

BlossomFinder::Blossom BlossomFinder::calculateBlossomFor(const Graph::NodeMap <size_t>& nodeToUF,
														  UnionFind& components,
														  size_t xIndex, double cutCost,
														  const Graph::NodeMap <size_t>& adjacentEDash) {
	//Die Blüte zur aktuellen Kante
	Blossom curr;
	//Informationen zur Kante e mit abs(c'(e)-c(e)) bzw. abs(c(e)-0.5) minimal:
	//Die Kante selbst
	Graph::Edge minDiffEdge = lemon::INVALID;
	//c(e)-0.5
	double minDiffVal = std::numeric_limits<double>::max();

	/*
	 * Der Index der Kante in curr.teeth, um sie schnell entfernen zu können, oder std::numeric_limits<size_t>::max(),
	 * falls die Kante nicht in curr.teeth enthalten ist.
	 */
	size_t minIndex = std::numeric_limits<size_t>::max();
	//Zähne ohne Beachtung der Parität berechnen
	for (Graph::EdgeIt eIt(fractionalGraph); eIt != lemon::INVALID; ++eIt) {
		bool uInX = components.find(nodeToUF[fractionalGraph.u(eIt)]) == xIndex;
		bool vInX = components.find(nodeToUF[fractionalGraph.v(eIt)]) == xIndex;
		//Ein Ende ist in der "unteren" Komponente, das andere nicht
		if (uInX != vInX) {
			bool inF = capacitiesFractional[eIt] > 0.5;
			if (inF) {
				curr.teeth.push_back(eIt);
			}
			//Neue minimale Kante gefunden
			if (std::abs(capacitiesFractional[eIt] - 0.5) < std::abs(minDiffVal)) {
				minDiffEdge = eIt;
				minDiffVal = capacitiesFractional[eIt] - 0.5;
				if (inF) {
					minIndex = curr.teeth.size() - 1;
				} else {
					minIndex = std::numeric_limits<size_t>::max();
				}
			}
		}
	}
	//Die Kardinalität von Xf geschnitten mit T'
	unsigned xAndTDash = 0;
	for (Graph::NodeIt nIt(fractionalGraph); nIt != lemon::INVALID; ++nIt) {
		//Ist der Knoten im Griff?
		if (components.find(nodeToUF[nIt]) == xIndex) {
			//Ist der Knoten in T delta V'?
			if (oddNodes[nIt] != (adjacentEDash[nIt] % 2 == 1)) {
				++xAndTDash;
			}
			curr.handle.push_back(nIt);
		}
	}
	//Ist der Schnitt gültig?
	bool valid = true;
	if (xAndTDash % 2 == 0) {
		if (minDiffEdge == lemon::INVALID) {
			//Die Paritätsbedingung kann nicht erfüllt werden, es gibt keine Kante im Schnitt
			valid = false;
		} else if (minIndex < curr.teeth.size()) {
			//Kante aus den Zähnen entfernen
			cutCost += 2 * minDiffVal;
			curr.teeth[minIndex] = curr.teeth.back();
			curr.teeth.pop_back();
		} else {
			//Kante zu den Zähnen hinzufügen
			cutCost -= 2 * minDiffVal;
			curr.teeth.push_back(minDiffEdge);
		}
	}
	//Falls die Blüte gültig ist und Wert kleiner als 1 hat, wird sie zurückgegeben
	if (valid && tolerance.less(cutCost, 1)) {
		return curr;
	} else {
		return {};
	}
}

/**
 * Erzeugt den fraktionalen Graphen zum Eingabegraphen
 */
void BlossomFinder::setupFractionalGraph() {
	//Ordnet jedem Knoten den entsprechenden Knoten im fraktionalen Graphen zu
	Graph::NodeMap <Graph::Node> mainToFract(inputGraph, lemon::INVALID);
	for (Graph::NodeIt it(inputGraph); it != lemon::INVALID; ++it) {
		Graph::Node newNode = fractionalGraph.addNode();
		contraction[newNode] = {it};
		mainToFract[it] = newNode;
	}
	for (Graph::EdgeIt it(inputGraph); it != lemon::INVALID; ++it) {
		Graph::Node u = inputGraph.u(it);
		Graph::Node v = inputGraph.v(it);
		if (!tolerance.different(capacitiesMain[it], 1)) {
			//1-Kante entfernen
			oneEdges.push_back(it);
			oddNodes[mainToFract[u]] = !oddNodes[mainToFract[u]];
			oddNodes[mainToFract[v]] = !oddNodes[mainToFract[v]];
		} else {
			Graph::Edge newEdge = fractionalGraph.addEdge(mainToFract[u], mainToFract[v]);
			fractToMainEdges[newEdge] = it;
			capacitiesFractional[newEdge] = capacitiesMain[it];
		}
	}
}

/**
 * Ersetzt eine Blüte in einem durch Kontraktion entstandenen Graphen durch die Blüte im ursprünglichen Graphen
 * @param nodeMap Ordnet jedem Knoten die entsprechenden ursprünglichen Knoten zu
 * @param edgeMap Ordnet jeder Kante die entsprechende ursprüngliche Kante zu
 */
void BlossomFinder::Blossom::replaceByMapped(const tsp_util::ContractionMapGraph& nodeMap,
											 const Graph::EdgeMap <Graph::Edge>& edgeMap) {
	std::vector<Graph::Node> newHandle;
	for (Graph::Node& n:handle) {
		newHandle.insert(newHandle.end(), nodeMap[n].begin(), nodeMap[n].end());
	}
	handle = newHandle;
	for (Graph::Edge& e:teeth) {
		e = edgeMap[e];
	}
}

/**
 * @return Ob die Blüte eine "echte" Blüte mit mindestens 3 Zähnen ist, die eine 2-Matching-Constraint erzeugt, oder
 * eine "falsche" Blüte mit nur einem Zahn, die eine Subtour-Constraint erzeugt
 */
bool BlossomFinder::Blossom::isProperBlossom() const {
	return teeth.size() >= 3;
}
