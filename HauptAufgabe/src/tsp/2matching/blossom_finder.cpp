#include <blossom_finder.hpp>
#include <lemon/gomory_hu.h>

BlossomFinder::BlossomFinder(const Graph& g, Graph::EdgeMap<double>& capacities, lemon::Tolerance<double> tolerance,
							 bool contractPaths) :
		mainGraph(g), tolerance(tolerance), capacitiesMain(capacities), capacitiesFractional(fractionalGraph),
		shouldContractPaths(contractPaths), nodeCount(lemon::countNodes(mainGraph)) {}

std::vector<BlossomFinder::Blossom> BlossomFinder::findViolatedBlossoms() {
	fractionalGraph.clear();
	tsp_util::ContractionMapGraph contraction(fractionalGraph);
	Graph::EdgeMap <Graph::Edge> fractToMainEdges(fractionalGraph, lemon::INVALID);
	Graph::NodeMap<bool> odd(fractionalGraph, false);
	std::vector<Graph::Edge> oneEdges;
	{
		Graph::NodeMap <Graph::Node> mainToFract(mainGraph, lemon::INVALID);
		for (Graph::NodeIt it(mainGraph); it != lemon::INVALID; ++it) {
			Graph::Node newNode = fractionalGraph.addNode();
			contraction[newNode] = {it};
			mainToFract[it] = newNode;
			for (Graph::OutArcIt arcIt(mainGraph, it); arcIt != lemon::INVALID; ++arcIt) {
				Graph::Node target = mainGraph.target(arcIt);
				if (mainToFract[target] != lemon::INVALID) {
					if (!tolerance.different(capacitiesMain[arcIt], 1)) {
						oneEdges.push_back(arcIt);
						odd[newNode] = !odd[newNode];
						odd[mainToFract[target]] = !odd[mainToFract[target]];
					} else {
						Graph::Edge newEdge = fractionalGraph.addEdge(newNode, mainToFract[target]);
						fractToMainEdges[newEdge] = arcIt;
						capacitiesFractional[newEdge] = capacitiesMain[arcIt];
					}
				}
			}
		}
	}
	if (shouldContractPaths) {
		contractPaths(odd, contraction);
	}
	std::vector<Blossom> violated = lemma1220(odd);
	for (size_t i = 0; i < violated.size();) {
		Blossom& b = violated[i];
		b.replaceByMapped(contraction, fractToMainEdges);
		if (finalizeBlossom(b, oneEdges)) {
			++i;
			assert(b.teeth.size() % 2 == 1);
		} else {
			std::swap(b, violated.back());
			violated.pop_back();
		}
	}
	return violated;
}

/**
 * Berechnet wie in Lemma 12.20 Blüten mit Wert kleiner als 1. c'(e) ist für alle Kanten 1-c(e)
 * @param graph Der zu betrachtende Graph
 * @param odd die ungeraden Knoten, bzw. die Menge T
 * @param c Die Kostenfunktion c
 * @return ein Vector mit Blüten mit Wert kleiner als 1. Falls es solche gibt, ist der Vector nicht leer.
 */
std::vector<BlossomFinder::Blossom> BlossomFinder::lemma1220(const Graph::NodeMap<bool>& odd) {
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
			Blossom b = calculateAndAddBlossom(nodeToUF, components, xIndex, cutCost, odd, adjacentEDash);
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
void BlossomFinder::contractPaths(Graph::NodeMap<bool>& odd, tsp_util::ContractionMapGraph& toOrig) {
	bool contracted;
	do {
		contracted = false;
		for (Graph::NodeIt start(fractionalGraph); start != lemon::INVALID && !contracted; ++start) {
			if (odd[start]) {
				contracted = findAndContractPath(start, toOrig, odd);
			}
		}
	} while (contracted);
}

/**
 * Findet und kontrahiert einen start enthaltenden alternierenden Pfad in G, falls es einen solchen gibt
 * @return true, genau dann wenn ein Pfad kontrahiert wurde
 */
bool BlossomFinder::findAndContractPath(Graph::Node start, tsp_util::ContractionMapGraph& toOrig,
										Graph::NodeMap<bool>& odd) {
	if (!isValidInternalNode(fractionalGraph, start, capacitiesFractional)) {
		return false;
	}
	Graph::NodeMap<bool> inPath(fractionalGraph);
	//Wird im ersten discoverPath-Aufruf auf die erste Kante des linken Teilpfades gesetzt
	Graph::Edge firstEdge = lemon::INVALID;
	//Teilpfad in eine Richtung ("links") finden
	std::vector<Graph::Node> left = discoverPath(start, firstEdge, odd, inPath);
	std::vector<Graph::Node> right;
	//Ein isolierter Kreis aus ungeraden Knoten
	bool isCycle = left.front() == left.back();
	//Falls der Pfad noch kein Kreis ist, gibt es einen rechten Teilpfad
	if (!isCycle) {
		//Anderen Teilpfad finden
		right = discoverPath(start, firstEdge, odd, inPath);
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
	bool resultOdd = odd[remainingNode];
	std::vector<Graph::Node>& contractedSet = toOrig[remainingNode];
	//Kontrahieren
	for (Graph::Node remove:toRemove) {
		assert(remove != remainingNode);
		if (odd[remove]) {
			resultOdd = !resultOdd;
		}
		contractedSet.insert(contractedSet.end(), toOrig[remove].begin(), toOrig[remove].end());
		fractionalGraph.erase(remove);
	}
	odd[remainingNode] = resultOdd;
	return true;
}

/**
 * Gibt den Wert einer Kante in einem alternierenden Pfad an, der start als inneren Knoten enthält, oder -1, falls es
 * keinen solchen gibt
 * @param g Der Graph, in dem der alternierende Pfad existieren soll
 * @param start der Knoten, der im Pfad enthalten sein soll
 * @param c Die Kantengewichte
 * @return das gewicht einer Kante im alternierenden Pfad oder -1, falls es keinen gibt
 */
bool BlossomFinder::isValidInternalNode(const Graph& g, Graph::Node start, const Graph::EdgeMap<double>& c) {
	size_t degree = 0;
	double edgeSum = 0;
	for (Graph::OutArcIt it(g, start); it != lemon::INVALID && degree <= 2 && !tolerance.less(1, edgeSum); ++it) {
		++degree;
		edgeSum += c[it];
	}
	return degree == 2 && !tolerance.different(edgeSum, 1);
}

/**
 * Findet einen Teil eines alternierenden Pfades
 * @param graph der Graph, in dem der alternierende Pfad gefunden werden soll
 * @param start der Ausgangsknoten. Muss als innerer Knoten eines alternierenden Pfades gültig sein.
 * @param exclude Entweder lemon::INVALID: in diesem Fall wird die erste Kante im Teilpfad hier gespeichert; oder eine
 * Kante, die nicht im Teilpfad genutzt werden darf
 * @param odd Die ungeraden Knoten
 * @param visited Gibt an, ob ein Knoten schon im Pfad enthalten ist
 * @param c Die Kantengewichte
 * @return den Teilpfad eines alternierenden Pfades mit den angegebenen Eigenschaften. Das erste Element ist immer start
 */
std::vector<Graph::Node> BlossomFinder::discoverPath(Graph::Node start, lemon::ListGraphBase::Edge& exclude,
													 const Graph::NodeMap<bool>& odd, Graph::NodeMap<bool>& visited) {
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
			//TODO warum gibt it!=lastEdge eigentlich einen Compiler-Error?
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
		canContinuePath = fractionalGraph.valid(nextNode) && !visited[nextNode] && odd[nextNode]
						  && isValidInternalNode(fractionalGraph, nextNode, capacitiesFractional);
		current = nextNode;
	} while (canContinuePath);
	if (current != lemon::INVALID) {
		path.push_back(current);
		visited[current] = true;
	}
	return path;
}

//TODO better name?
bool BlossomFinder::finalizeBlossom(Blossom& b, const std::vector<Graph::Edge>& oneEdges) {
	std::vector<Graph::Edge>& teeth = b.teeth;
	std::vector<Graph::Node>& handle = b.handle;
	if (handle.size() < 2 || handle.size() > nodeCount - 2) {
		return false;
	}
	Graph::NodeMap<bool> inHandle(mainGraph, false);
	for (Graph::Node n:handle) {
		inHandle[n] = true;
	}
	//Gibt an, ob eine Variable/Kante im TSP-Graphen in F ist
	Graph::EdgeMap<bool> fTSP(mainGraph, false);
	//Die Menge F vor dem Umwandeln zu einem Matching
	//Alle 1-Kanten im Schnitt von X sind in F
	for (Graph::Edge e:oneEdges) {
		if (inHandle[mainGraph.u(e)] != inHandle[mainGraph.v(e)]) {
			teeth.push_back(e);
		}
	}
	assert(teeth.size() % 2 == 1);
	//Die Kanten im berechneten F zu F hinzufügen
	for (Graph::Edge tooth:teeth) {
		fTSP[tooth] = true;
	}
	{
		Graph::NodeMap <size_t> incident(mainGraph);
		for (Graph::Edge tooth:teeth) {
			++incident[mainGraph.u(tooth)];
			if (incident[mainGraph.u(tooth)] >= 3) {
				return false;
			}
			++incident[mainGraph.v(tooth)];
			if (incident[mainGraph.v(tooth)] >= 3) {
				return false;
			}
		}
	}
	size_t handleSize = b.handle.size();
	//Ordnet jedem Knoten die inzidente Kante in F zu
	Graph::NodeMap <Graph::Edge> incidentF(mainGraph, lemon::INVALID);
	for (Graph::Edge e:teeth) {
		Graph::Node ends[2] = {mainGraph.u(e), mainGraph.v(e)};
		for (Graph::Node end:ends) {
			if (incidentF[end] == lemon::INVALID) {
				//Es ist noch keine Kante in F zu end inzident
				incidentF[end] = e;
			} else {
				//Die Kanten entfernen, die ein gemeinsames Ende haben
				Graph::Edge oldEdge = incidentF[end];
				incidentF[mainGraph.u(oldEdge)] = lemon::INVALID;
				incidentF[mainGraph.v(oldEdge)] = lemon::INVALID;
				incidentF[mainGraph.u(e)] = lemon::INVALID;
				incidentF[mainGraph.v(e)] = lemon::INVALID;
				fTSP[oldEdge] = false;
				fTSP[e] = false;
				//Das gemeinsame Ende aus X entfernen bzw zu X hinzufügen
				inHandle[end] = !inHandle[end];
				if (inHandle[end]) {
					++handleSize;
				} else {
					--handleSize;
				}
				break;
			}
		}
	}
	size_t i = 0;
	while (i < teeth.size()) {
		Graph::Edge tooth = teeth[i];
		if (!fTSP[tooth]) {
			std::swap(teeth[i], teeth.back());
			teeth.pop_back();
		} else {
			assert(inHandle[mainGraph.u(tooth)] != inHandle[mainGraph.v(tooth)]);
			++i;
		}
	}
	/*
	 * Wahr, falls das berechnete X verwendet werden soll; falsch, falls das Komplement verwendet werden soll,
	 * um eine dünner besetzte Constraint zu erhalten.
	 */
	//TODO
	//const bool handleVal = handleSize < nodeCount / 2;
	//handle.clear();
	//for (Graph::NodeIt it(mainGraph); it != lemon::INVALID; ++it) {
	//	if (inHandle[it] == handleVal) {
	//		handle.push_back(it);
	//	}
	//}
	return handle.size() >= 2 && handle.size() <= nodeCount - 2;
}

BlossomFinder::Blossom BlossomFinder::calculateAndAddBlossom(const Graph::NodeMap <size_t>& nodeToUF,
															 UnionFind& components, size_t xIndex,
															 double cutCost, const Graph::NodeMap<bool>& odd,
															 const Graph::NodeMap <size_t>& adjacentEDash) {
	//Die Blüte zur aktuellen Kante
	Blossom curr;
	//Informationen zur Kante mit abs(c'(e)-c(e)) bzw. abs(c(e)-0.5) minimal:
	//Die Kante selbst
	Graph::Edge minDiffEdge = lemon::INVALID;
	//c(e)-0.5
	double minDiffVal = std::numeric_limits<double>::max();
	/*
	 * Der Index der Kante in curr.teeth, um sie schnell entfernen zu können, oder std::numeric_limits<size_t>::max(),
	 * falls die Kante nicht in curr.teeth enthalten ist.
	 */
	size_t minIndex = std::numeric_limits<size_t>::max();
	//F ohne Beachtung der Parität berechnen
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
		//Ist der Knoten in X?
		if (components.find(nodeToUF[nIt]) == xIndex) {
			//Ist der Knoten in T delta V'?
			if (odd[nIt] != (adjacentEDash[nIt] % 2 == 1)) {
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
			//Kante aus F entfernen
			cutCost += 2 * minDiffVal;
			curr.teeth[minIndex] = curr.teeth.back();
			curr.teeth.pop_back();
		} else {
			//Kante zu F hinzufügen
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

bool BlossomFinder::Blossom::isProperBlossom() const {
	return teeth.size() >= 3;
}
