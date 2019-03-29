#include <two_matching_cut_gen.hpp>
#include <lemon/unionfind.h>
#include <union_find.hpp>
#include <tsp_utils.hpp>

TwoMatchingCutGen::TwoMatchingCutGen(const TSPInstance& inst, bool contract)
		: tsp(inst), enableContraction(contract), tolerance(1e-5) {}

/**
 * Findet verletzte 2-Matching-Constraints, falls es solche gibt. Der Algorithmus entspricht grob dem aus Satz 12.21 mit
 * einigen Änderungen:
 * 1. z ist hier nicht nötig: Die Grad-Constraints sind immer erfüllt, also gilt für alle an z anliegenden Kanten e
 * c(e)=0 und c'(e)=inf. Sie können also nie in der Menge F einer verletzten 2-Matching-Constraint enthalten sein.
 * 2. Der Graph wird vor dem Anwenden von Lemma 12.20 durch Pfad-Kontraktion nach Proposition 4.6 in
 * Padberg, M. & Rinaldi, G. Mathematical Programming (1990) 47: 219. https://doi.org/10.1007/BF01580861
 * 3. Im Algorithmus aus Lemma 12.21 werden hier nicht die minimalen Blüten berechnet, sondern alle aus dem Gomory-Hu-
 * Baum entstehenden Blüten mit Wert echt kleiner 1.
 * 4. Die Blüten werden vor dem Hinzufügen der Constraints so verändert, dass F ein Matching ist. Der Algorithmus wird
 * im selben Paper wie oben ohne Beweis gegeben. TODO: Quelle mit Beweis finden oder selbst einen schreiben
 */
CutGenerator::CutStatus TwoMatchingCutGen::validate(LinearProgram& lp, const std::vector<double>& solution,
													CutStatus currentStatus) {
	if (currentStatus==CutGenerator::maybe_recalc) return CutGenerator::valid;
	Graph workGraph;
	//Ordnet den Städten der TSP-Instanz einen Knoten im Arbeitsgraphen zu. Nur bis zum Aufruf von contractPaths gültig!
	std::vector<Graph::Node> origToWork(static_cast<size_t>(tsp.getCityCount()));
	//Ordnet den Knoten im Arbeitsgraphen die Städte in der TSP-Instanz zu. Auch nach dem Aufruf von contractPaths gültig
	ContractionMap workToOrig(workGraph);
	//Ungerade Knoten, bzw. Knoten in T
	Graph::NodeMap<bool> odd(workGraph);
	//Ordnet den Kanten die zugehörigen Variablen zu
	Graph::EdgeMap <variable_id> toVariable(workGraph);
	//c wie in Satz 12.21
	Graph::EdgeMap<double> c(workGraph);
	Graph::NodeMap <city_id> workToOrigTemp(workGraph);
	//Die Kanten mit Wert 1
	std::vector<variable_id> oneEdges = tsp_util::createFractionalGraph(tsp, tolerance, solution, workGraph,
																		workToOrigTemp, origToWork, toVariable, c, odd);
	for (Graph::NodeIt it(workGraph); it!=lemon::INVALID; ++it) {
		workToOrig[it] = {workToOrigTemp[it]};
	}
	if (enableContraction) {
		contractPaths(workGraph, odd, c, workToOrig);
	}
	std::vector<Blossom> allMin = lemma1220(workGraph, odd, c);
	if (!allMin.empty()) {
		std::vector<LinearProgram::Constraint> constrs;
		constrs.reserve(allMin.size());
		for (Blossom& min:allMin) {
			//Gibt an, ob ein Knoten in der aktuellen Menge X ist
			std::vector<bool> isInX(static_cast<size_t>(tsp.getCityCount()));
			//Die Größe der Menge X
			size_t sizeX = 0;
			for (Graph::Node contracted:min.x) {
				for (city_id orig:workToOrig[contracted]) {
					isInX[orig] = true;
				}
				sizeX += workToOrig[contracted].size();
			}
			std::vector<variable_id> f;
			f.reserve(min.f.size());
			for (Graph::Edge e:min.f) {
				f.push_back(toVariable[e]);
			}
			finalizeBlossom(oneEdges, isInX, f, sizeX);
			const size_t sizeF = f.size();
			/*
			 * Mit |F|==1 ist auch die Subtour-Constraint für X verletzt und impliziert die
			 * 2-Matching-Constraint
			 */
			assert(sizeF%2==1);
			//Die Indizes der Variablen in den hinzugefügten Constraints
			std::vector<variable_id> indices;
			if (sizeF>1) {
				for (variable_id e:f) {
					indices.push_back(e);
					assert(isInX[tsp.getLowerEnd(e)]!=isInX[tsp.getHigherEnd(e)]);
				}
			}
			std::vector<city_id> xElements;
			/*
			 * Wahr, falls das berechnete X verwendet werden soll; falsch, falls das Komplement verwendet werden soll,
			 * um eine dünner besetzte Constraint zu erhalten.
			 */
			const bool valForX = sizeX<tsp.getCityCount()/2;
			//Die Kanten des induzierten Graphen hinzufügen
			for (city_id i = 0; i<tsp.getCityCount(); ++i) {
				if (isInX[i]==valForX) {
					for (city_id other:xElements) {
						indices.push_back(tsp.getVariable(i, other));
					}
					xElements.push_back(i);
				}
			}
			if (!indices.empty()) {
				if (sizeF > 1) {
					//2-Matching-Constraint
					constrs.emplace_back(indices, std::vector<double>(indices.size(), 1), LinearProgram::less_eq,
										 static_cast<size_t>(xElements.size() + sizeF / 2));
				} else {
					//Subtour-Constraint
					constrs.emplace_back(indices, std::vector<double>(indices.size(), 1), LinearProgram::less_eq,
										 xElements.size() - 1);
				}
			}
		}
		if (!constrs.empty()) {
			lp.addConstraints(constrs);
			return CutGenerator::maybe_recalc;
		}
	}
	//if (enableContraction) {
	//	TwoMatchingCutGen tmp(tsp, false);
	//	assert(tmp.validate(lp, solution)==CutGenerator::valid);
	//}
	return CutGenerator::valid;
}

/**
 * Berechnet wie in Lemma 12.20 Blüten mit Wert kleiner als 1. c'(e) ist für alle Kanten 1-c(e)
 * @param graph Der zu betrachtende Graph
 * @param odd die ungeraden Knoten, bzw. die Menge T
 * @param c Die Kostenfunktion c
 * @return ein Vector mit Blüten mit Wert kleiner als 1. Falls es solche gibt, ist der Vector nicht leer.
 */
std::vector<TwoMatchingCutGen::Blossom> TwoMatchingCutGen::lemma1220(const Graph& graph,
																	 const Graph::NodeMap<bool>& odd,
																	 const Graph::EdgeMap<double>& c) {
	//d wie im Beweis von Lemma 12.20
	Graph::EdgeMap<double> d(graph);
	//Die Anzahl der inzidenten Kanten in E'
	Graph::NodeMap <size_t> adjacentEDash(graph);
	//d und adjacentEDash berechnen
	for (Graph::EdgeIt it(graph); it!=lemon::INVALID; ++it) {
		if (c[it]>0.5) {
			++adjacentEDash[graph.v(it)];
			++adjacentEDash[graph.u(it)];
			d[it] = 1-c[it];
		} else {
			d[it] = c[it];
		}
	}
	//Gomory-Hu-Baum/Arboreszenz bezüglich d berechnen
	lemon::GomoryHu<Graph, Graph::EdgeMap<double>>
	gh(graph, d);
	gh.run();
	/*
	 * Die Arboreszenz durchlaufen. Eine einfachere Implementierung, bei der wiederholt minCutMap aufgerufen wird, ist
	 * deutlich langsamer. Die Arboreszenz wird von den Blättern aus zur Wurzel durchlaufen. Die Knoten auf der "unteren"
	 * Seite der aktuellen Kante sind als Menge in components gespeichert, nach dem Bearbeiten einer Kante werden die
	 * Komponenten an den Enden der Kante zusammengefügt.
	 */
	//Die Anzahl der Kinder eines Knotens in der Gomory-Hu-Arboreszenz
	Graph::NodeMap <size_t> childCount(graph);
	//Die Anzahl der Knoten im Graphen bzw. in der Arboreszenz
	size_t nodeCount = 0;
	for (Graph::NodeIt it(graph); it!=lemon::INVALID; ++it) {
		Graph::Node pred = gh.predNode(it);
		if (pred!=lemon::INVALID) {
			++childCount[pred];
		}
		++nodeCount;
	}
	//Ordnet jedem Knoten eine ID im Union-Find zu
	Graph::NodeMap <size_t> nodeToUF(graph);
	UnionFind components(nodeCount);
	/*
	 * Die Blätter der Arboreszenz. Wenn eine Kante behandelt wurde, wird ihr "unterer" Endknoten aus der Arboreszenz
	 * "entfernt", der "obere" Endknoten wird also evtl zu einem Blatt und wird zu diesem Vector hinzugefügt.
	 */
	std::vector<Graph::Node> leaves;
	{
		size_t nextIndex = 0;
		for (Graph::NodeIt it(graph); it!=lemon::INVALID; ++it, ++nextIndex) {
			if (childCount[it]==0) {
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
		if (pred==lemon::INVALID) {
			continue;
		}
		//Repräsentant der "unteren" Seite der aktuellen Kante
		const size_t xIndex = components.find(nodeToUF[currentNode]);
		double cutCost = gh.predValue(currentNode);
		//Der Wert der Blüte ist mindestens der Wert des Schnitts
		if (tolerance.less(cutCost, 1)) {
			Blossom b = calculateAndAddBlossom(nodeToUF, components, xIndex, graph, cutCost, c, odd, adjacentEDash);
			if (!b.x.empty()) {
				ret.push_back(b);
			}
		}
		components.mergeRoots(xIndex, components.find(nodeToUF[pred]));
		--childCount[pred];
		if (childCount[pred]==0) {
			leaves.push_back(pred);
		}
	}
	return ret;
}

/**
 * Kontrahiert so lange alternierende Pfade in G, bis es keine mehr gibt, und speichert die Kontraktionen in toOrig
 */
void TwoMatchingCutGen::contractPaths(Graph& g, Graph::NodeMap<bool>& odd, Graph::EdgeMap<double>& c,
									  ContractionMap& toOrig) {
	bool contracted;
	do {
		contracted = false;
		for (Graph::NodeIt start(g); start!=lemon::INVALID && !contracted; ++start) {
			if (odd[start]) {
				contracted = findAndContractPath(g, start, toOrig, odd, c);
			}
		}
	} while (contracted);
}

/**
 * Findet und kontrahiert einen start enthaltenden alternierenden Pfad in G, falls es einen solchen gibt
 * @return true, genau dann wenn ein Pfad kontrahiert wurde
 */
bool TwoMatchingCutGen::findAndContractPath(Graph& g, Graph::Node start, ContractionMap& toOrig,
											Graph::NodeMap<bool>& odd, const Graph::EdgeMap<double>& c) {
	if (!isValidInternalNode(g, start, c)) {
		return false;
	}
	Graph::NodeMap<bool> inPath(g);
	//Wird im ersten discoverPath-Aufruf auf die erste Kante des linken Teilpfades gesetzt
	Graph::Edge firstEdge = lemon::INVALID;
	//Teilpfad in eine Richtung ("links") finden
	std::vector<Graph::Node> left = discoverPath(g, start, firstEdge, odd, inPath, c);
	std::vector<Graph::Node> right;
	//Ein isolierter Kreis aus ungeraden Knoten
	bool isCycle = left.front()==left.back();
	//Falls der Pfad noch kein Kreis ist, gibt es einen rechten Teilpfad
	if (!isCycle) {
		//Anderen Teilpfad finden
		right = discoverPath(g, start, firstEdge, odd, inPath, c);
		//Kreis mit einem geraden Knoten
		isCycle = right.back()==left.back();
	}
	//Der Knoten, der von der kontrahierten Menge übrig bleibt
	Graph::Node remainingNode = left.back();
	//Die restlichen Knoten der kontrahierten Menge:
	//Der letzte Knoten von left wird nie kontrahiert (er bleibt übrig), der erste Knoten ist entweder auch der letzte
	//Knoten oder gleich dem ersten Knoten der rechten Seite
	std::vector<Graph::Node> toRemove(left.begin()+1, left.end()-1);
	if (right.size()>1) {
		//Der letzte rechte Knoten wird nicht hinzugefügt: Falls der Pfad ein Kreis ist, ist dies der verbleibende Knoten;
		//sonst kann nur P-u_1 kontrahiert werden
		toRemove.insert(toRemove.end(), right.begin(), right.end()-1);
	}
	if (!isCycle) {
		//Falls der Pfad kein Kreis ist, muss die Kante zwischen dem letzten kontrahierten Knoten und u_1 verbleiben
		//Es kann nicht einfach die Kante gelöscht und durch eine neue ersetzt werden, da sonst Werte in EdgeMap's
		//verloren gehen (z.B. toVariable in validate)
		Graph::Node lastContracted = right[right.size()-2];
		Graph::Node otherRemaining = right.back();
		for (Graph::OutArcIt it(g, lastContracted); it!=lemon::INVALID; ++it) {
			Graph::Node target = g.target(it);
			if (target==otherRemaining) {
				if (g.u(it)==lastContracted) {
					g.changeU(it, remainingNode);
				} else {
					g.changeV(it, remainingNode);
				}
				break;
			}
		}
	}
	//Ob der durch die Kontraktion entstehende Knoten ungerade ist
	bool resultOdd = odd[remainingNode];
	std::vector<city_id>& contractedSet = toOrig[remainingNode];
	//Kontrahieren
	for (Graph::Node remove:toRemove) {
		assert(remove!=remainingNode);
		if (odd[remove]) {
			resultOdd = !resultOdd;
		}
		contractedSet.insert(contractedSet.end(), toOrig[remove].begin(), toOrig[remove].end());
		g.erase(remove);
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
bool TwoMatchingCutGen::isValidInternalNode(const Graph& g, Graph::Node start, const Graph::EdgeMap<double>& c) {
	size_t degree = 0;
	double edgeSum = 0;
	for (Graph::OutArcIt it(g, start); it!=lemon::INVALID && degree<=2 && !tolerance.less(1, edgeSum); ++it) {
		++degree;
		edgeSum += c[it];
	}
	return degree==2 && !tolerance.different(edgeSum, 1);
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
std::vector<Graph::Node> TwoMatchingCutGen::discoverPath(const Graph& graph, Graph::Node start,
														 lemon::ListGraphBase::Edge& exclude,
														 const Graph::NodeMap<bool>& odd, Graph::NodeMap<bool>& visited,
														 const Graph::EdgeMap<double>& c) {
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
		for (Graph::OutArcIt it(graph, current); it!=lemon::INVALID && nextNode==lemon::INVALID; ++it) {
			//TODO warum gibt it!=lastEdge eigentlich einen Compiler-Error?
			//current ist ein gültiger innerer Knoten, hat also Grad 2 und Kantensumme 1. it ist nicht die vorherige
			//Kante im Pfad, also muss it die nächste Kante des Pfades sein
			if (lastEdge!=it) {
				nextNode = graph.target(it);
				lastEdge = it;
				//Erster Fall der Definition von exclude
				if (exclude==lemon::INVALID) {
					exclude = it;
				}
			}
		}
		//Prüfen, ob der nächste Knoten ein gültiger innerer Knoten ist
		assert(graph.valid(nextNode));
		canContinuePath = graph.valid(nextNode) && !visited[nextNode] && odd[nextNode]
						  && isValidInternalNode(graph, nextNode, c);
		current = nextNode;
	} while (canContinuePath);
	if (current!=lemon::INVALID) {
		path.push_back(current);
		visited[current] = true;
	}
	return path;
}

//TODO better name?
void TwoMatchingCutGen::finalizeBlossom(const std::vector<variable_id>& oneEdges, std::vector<bool>& isInX,
										std::vector<variable_id>& f, size_t& sizeX) {
//Gibt an, ob eine Variable/Kante im TSP-Graphen in F ist
	std::vector<bool> fTSP(static_cast<size_t>(tsp.getEdgeCount()));
	//Die Menge F vor dem Umwandeln zu einem Matching
	//Alle 1-Kanten im Schnitt von X sind in F
	for (variable_id e:oneEdges) {
		city_id endU = tsp.getLowerEnd(e);
		city_id endV = tsp.getHigherEnd(e);
		if (isInX[endU]!=isInX[endV]) {
			fTSP[e] = true;
			f.push_back(e);
		}
	}
	//Die Kanten im berechneten F zu F hinzufügen
	for (variable_id e:f) {
		fTSP[e] = true;
	}
	//Ordnet jedem Knoten die inzidente Kante in F zu
	std::vector<variable_id> incidentF(tsp.getCityCount(), LinearProgram::invalid_variable);
	for (variable_id e:f) {
		city_id ends[2] = {tsp.getLowerEnd(e), tsp.getHigherEnd(e)};
		for (city_id end:ends) {
			if (incidentF[end]==LinearProgram::invalid_variable) {
				//Es ist noch keine Kante in F zu end inzident
				incidentF[end] = e;
			} else {
				//Die Kanten entfernen, die ein gemeinsames Ende haben
				variable_id oldEdge = incidentF[end];
				incidentF[tsp.getLowerEnd(oldEdge)] = LinearProgram::invalid_variable;
				incidentF[tsp.getHigherEnd(oldEdge)] = LinearProgram::invalid_variable;
				incidentF[tsp.getLowerEnd(e)] = LinearProgram::invalid_variable;
				incidentF[tsp.getHigherEnd(e)] = LinearProgram::invalid_variable;
				fTSP[oldEdge] = false;
				fTSP[e] = false;
				//Das gemeinsame Ende aus X entfernen bzw zu X hinzufügen
				isInX[end] = !isInX[end];
				//Die Größen von X und F aktualisieren
				if (isInX[end]) {
					++sizeX;
				} else {
					--sizeX;
				}
				break;
			}
		}
	}
	size_t i = 0;
	while (i<f.size()) {
		variable_id fElement = f[i];
		if (!fTSP[fElement]) {
			std::swap(f[i], f.back());
			f.pop_back();
		} else {
			++i;
		}
	}
}

TwoMatchingCutGen::Blossom TwoMatchingCutGen::calculateAndAddBlossom(const Graph::NodeMap <size_t>& nodeToUF,
																	 UnionFind& components,
																	 size_t xIndex, const Graph& g, double cutCost,
																	 const Graph::EdgeMap<double>& c,
																	 const Graph::NodeMap<bool>& odd,
																	 const Graph::NodeMap <size_t>& adjacentEDash) {
//Die Blüte zur aktuellen Kante
	Blossom curr;
	//Informationen zur Kante mit abs(c'(e)-c(e)) bzw. abs(c(e)-0.5) minimal:
	//Die Kante selbst
	Graph::Edge minDiffEdge = lemon::INVALID;
	//c(e)-0.5
	double minDiffVal = std::numeric_limits<double>::max();
	/*
	 * Der Index der Kante in curr.f, um sie schnell entfernen zu können, oder std::numeric_limits<size_t>::max(),
	 * falls die Kante nicht in curr.f enthalten ist.
	 */
	size_t minIndex = std::numeric_limits<size_t>::max();
	//F ohne Beachtung der Parität berechnen
	for (Graph::EdgeIt eIt(g); eIt!=lemon::INVALID; ++eIt) {
		bool uInX = components.find(nodeToUF[g.u(eIt)])==xIndex;
		bool vInX = components.find(nodeToUF[g.v(eIt)])==xIndex;
		//Ein Ende ist in der "unteren" Komponente, das andere nicht
		if (uInX!=vInX) {
			bool inF = c[eIt]>0.5;
			if (inF) {
				curr.f.push_back(eIt);
			}
			//Neue minimale Kante gefunden
			if (std::abs(c[eIt]-0.5)<std::abs(minDiffVal)) {
				minDiffEdge = eIt;
				minDiffVal = c[eIt]-0.5;
				if (inF) {
					minIndex = curr.f.size()-1;
				} else {
					minIndex = std::numeric_limits<size_t>::max();
				}
			}
		}
	}
	//Die Kardinalität von Xf geschnitten mit T'
	unsigned xAndTDash = 0;
	for (Graph::NodeIt nIt(g); nIt!=lemon::INVALID; ++nIt) {
		//Ist der Knoten in X?
		if (components.find(nodeToUF[nIt])==xIndex) {
			//Ist der Knoten in T delta V'?
			if (odd[nIt]!=(adjacentEDash[nIt]%2==1)) {
				++xAndTDash;
			}
			curr.x.push_back(nIt);
		}
	}
	//Ist der Schnitt gültig?
	bool valid = true;
	if (xAndTDash%2==0) {
		if (minDiffEdge==lemon::INVALID) {
			//Die Paritätsbedingung kann nicht erfüllt werden, es gibt keine Kante im Schnitt
			valid = false;
		} else if (minIndex<curr.f.size()) {
			//Kante aus F entfernen
			cutCost += 2*minDiffVal;
			curr.f[minIndex] = curr.f.back();
			curr.f.pop_back();
		} else {
			//Kante zu F hinzufügen
			cutCost -= 2*minDiffVal;
			curr.f.push_back(minDiffEdge);
		}
	}
	//Falls die Blüte gültig ist und Wert kleiner als 1 hat, wird sie zurückgegeben
	if (valid && tolerance.less(cutCost, 1)) {
		return curr;
	} else {
		return {};
	}
}
