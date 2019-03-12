#include <two_matching_cut_gen.hpp>
#include <lemon/unionfind.h>
#include <union_find.hpp>

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
CutGenerator::CutStatus TwoMatchingCutGen::validate(LinearProgram& lp, const std::vector<double>& solution) {
	Graph workGraph;
	//Ordnet den Städten der TSP-Instanz einen Knoten im Arbeitsgraphen zu. Nur bis zum Aufruf von contractPaths gültig!
	std::vector<Graph::Node> origToWork(static_cast<size_t>(tsp.getCityCount()));
	//Ordnet den Knoten im Arbeitsgraphen die Städte in der TSP-Instanz zu. Auch nach dem Aufruf von contractPaths gültig
	ContractionMap workToOrig(workGraph);
	for (city_id i = 0; i<tsp.getCityCount(); ++i) {
		Graph::Node newNode = workGraph.addNode();
		origToWork[i] = newNode;
		workToOrig[newNode] = {i};
	}
	//Ungerade Knoten, bzw. Knoten in T
	Graph::NodeMap<bool> odd(workGraph);
	//Ordnet den Kanten die zugehörigen Variablen zu
	Graph::EdgeMap <variable_id> toVariable(workGraph);
	//c wie in Satz 12.21
	Graph::EdgeMap<double> c(workGraph);
	//Die Kanten mit Wert 1
	std::vector<variable_id> oneEdges;
	for (city_id lower = 0; lower<tsp.getCityCount()-1; ++lower) {
		for (city_id higher = lower+1; higher<tsp.getCityCount(); ++higher) {
			variable_id varId = tsp.getVariable(higher, lower);
			/*
			 * Kanten mit Wert 0 können nicht in F enthalten sein, da der Wert der Blüte dann schon mindestens 1 wäre
			 * (und damit die Constraint nicht verletzt ist)
			 */
			if (tolerance.positive(solution[varId])) {
				Graph::Node endU = origToWork[lower];
				Graph::Node endV = origToWork[higher];
				if (tolerance.less(solution[varId], 1)) {
					//Kante hinzufügen
					Graph::Edge eWork = workGraph.addEdge(endU, endV);
					toVariable[eWork] = varId;
					c[eWork] = solution[varId];
				} else {
					/*
					 * Kanten im Schnitt mit Wert 1 müssen in F enthalten sein, sonst ist der Wert mindestens 1. Durch
					 * das Hinzufügen bzw Entfernen der Enden aus T bleiben die möglichen Blüten gleich.
					 */
					odd[endU] = !odd[endU];
					odd[endV] = !odd[endV];
					oneEdges.push_back(varId);
				}
			}
		}
	}
	if (enableContraction) {
		contractPaths(workGraph, odd, c, workToOrig);
	}
	std::vector<Blossom> allMin = lemma1220(workGraph, odd, c);
	if (!allMin.empty()) {
		//Die Indizes der Variablen in den hinzugefügten Constraints
		std::vector<variable_id> indices;
		//Die Indizes in indices, an denen die Constraints anfangen
		std::vector<int> constrStarts;
		//Die rechten Seiten der Constraints
		std::vector<double> rhs;
		for (Blossom& min:allMin) {
			//Gibt an, ob ein Knoten in der aktuellen Menge X ist
			std::vector<bool> isInX(tsp.getCityCount());
			//Die Größe der Menge X
			size_t sizeX = 0;
			for (Graph::Node contracted:min.x) {
				for (city_id orig:workToOrig[contracted]) {
					isInX[orig] = true;
				}
				sizeX += workToOrig[contracted].size();
			}
			//Gibt an, ob eine Variable/Kante im TSP-Graphen in F ist
			std::vector<bool> fTSP(tsp.getEdgeCount());
			//Die Menge F vor dem Umwandeln zu einem Matching
			std::vector<variable_id> prelimF;
			//Alle 1-Kanten im Schnitt von X sind in F
			for (variable_id e:oneEdges) {
				city_id endU = tsp.getLowerEnd(e);
				city_id endV = tsp.getHigherEnd(e);
				if (isInX[endU]!=isInX[endV]) {
					fTSP[e] = true;
					prelimF.push_back(e);
				}
			}
			//Die Kanten im berechneten F zu F hinzufügen
			for (Graph::Edge e:min.f) {
				fTSP[toVariable[e]] = true;
				prelimF.push_back(toVariable[e]);
			}
			//Ordnet jedem Knoten die inzidente Kante in F zu
			std::vector<variable_id> incidentF(tsp.getCityCount(), LinearProgram::invalid_variable);
			size_t sizeF = prelimF.size();
			for (variable_id e:prelimF) {
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
						sizeF -= 2;
						break;
					}
				}
			}
			constrStarts.push_back(static_cast<int>(indices.size()));
			/*
			 * Mit |F|==1 ist auch die Subtour-Constraint für X verletzt und impliziert die
			 * 2-Matching-Constraint
			 */
			if (sizeF>1) {
				for (variable_id e:prelimF) {
					if (fTSP[e]) {
						indices.push_back(e);
					}
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
			if (sizeF>1) {
				//2-Matching-Constraint
				rhs.push_back(static_cast<size_t>(xElements.size()+sizeF/2));
			} else {
				//Subtour-Constraint
				rhs.push_back(xElements.size()-1);
			}
		}
		lp.addConstraints(indices, std::vector<double>(indices.size(), 1), rhs, constrStarts,
						  std::vector<LinearProgram::CompType>(rhs.size(), LinearProgram::less_eq));
		return CutGenerator::maybe_recalc;
	} else {
		return CutGenerator::valid;
	}
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
	Graph::NodeMap <size_t> ufMap(graph);
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
			ufMap[it] = nextIndex;
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
		const size_t xIndex = components.find(ufMap[currentNode]);
		double cutCost = gh.predValue(currentNode);
		//Der Wert der Blüte ist mindestens der Wert des Schnitts
		if (tolerance.less(cutCost, 1)) {
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
			for (Graph::EdgeIt eIt(graph); eIt!=lemon::INVALID; ++eIt) {
				bool uInX = components.find(ufMap[graph.u(eIt)])==xIndex;
				bool vInX = components.find(ufMap[graph.v(eIt)])==xIndex;
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
			for (Graph::NodeIt nIt(graph); nIt!=lemon::INVALID; ++nIt) {
				//Ist der Knoten in X?
				if (components.find(ufMap[nIt])==xIndex) {
					//Ist der Knoten in T delta V'?
					if (odd[nIt]!=(adjacentEDash[nIt]%2==1)) {
						++xAndTDash;
					}
					curr.x.push_back(nIt);
				}
			}
			//Ist der Schnitt gültig?
			bool valid;
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
			} else {
				//Die Paritätsbedingung ist schon erfüllt
				valid = true;
			}
			//Falls die Blüte gültig ist und Wert kleiner als 1 hat, wird sie zurückgegeben
			if (valid && tolerance.less(cutCost, 1)) {
				ret.push_back(curr);
			}
		}
		components.mergeRoots(xIndex, components.find(ufMap[pred]));
		--childCount[pred];
		if (childCount[pred]==0) {
			leaves.push_back(pred);
		}
	}
	return ret;
}

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

bool TwoMatchingCutGen::findAndContractPath(Graph& g, Graph::Node start, ContractionMap& toOrig,
											Graph::NodeMap<bool>& odd,
											const Graph::EdgeMap<double>& c) {
	Graph::NodeMap<bool> inPath(g);
	double pathValue = -1;
	//Teilpfad in eine Richtung ("links") finden
	std::vector<Graph::Node> left = discoverPath(g, start, start, odd, inPath, c, pathValue);
	if (left.empty()) {
		return false;
	}
	pathValue = 1-pathValue;
	//Anderen Teilpfad finden
	std::vector<Graph::Node> right = discoverPath(g, start, left[0], odd, inPath, c, pathValue);
	//TODO when can right be empty?
	if (left.size()+right.size()<2 || right.empty()) {
		return false;
	}
	std::vector<Graph::Node> toContract(left.begin(), left.end()-1);
	if (!right.empty()) {
		toContract.insert(toContract.end(), right.begin(), right.end()-1);
	}
	Graph::Node remainingNode = left.back();
	if (remainingNode!=start) {
		toContract.push_back(start);
	}
	if (right.back()!=left.back()) {
		inPath[right.back()] = false;
		for (Graph::Node pathNode:toContract) {
			std::vector<Graph::Edge> toMove;
			for (Graph::OutArcIt out(g, pathNode); out!=lemon::INVALID; ++out) {
				Graph::Node target = g.target(out);
				if (!inPath[target]) {
					toMove.push_back(out);
				}
			}
			for (Graph::Edge e:toMove) {
				if (g.u(e)==pathNode) {
					g.changeU(e, remainingNode);
				} else {
					g.changeV(e, remainingNode);
				}
			}
		}
	}
	bool resultOdd = odd[remainingNode];
	std::vector<city_id>& contractedSet = toOrig[remainingNode];
	for (Graph::Node remove:toContract) {
		if (odd[remove]) {
			resultOdd = !resultOdd;
		}
		contractedSet.insert(contractedSet.end(), toOrig[remove].begin(), toOrig[remove].end());
		g.erase(remove);
	}
	odd[remainingNode] = resultOdd;
	return true;
}

std::vector<Graph::Node> TwoMatchingCutGen::discoverPath(const Graph& graph, Graph::Node start, Graph::Node exclude,
														 const Graph::NodeMap<bool>& odd, Graph::NodeMap<bool>& visited,
														 const Graph::EdgeMap<double>& c, double& firstEdgeVal) {
	//Die Knoten auf dem Pfad
	std::vector<Graph::Node> path;
	Graph::Node current = start;
	//Der Grad des aktuellen Knotens
	size_t degree;
	//Die Summe der Werte der inzidenten Kanten
	double edgeSum;
	//Der Wert der nächsten Kante
	double nextEdgeVal = firstEdgeVal;
	do {
		if (current!=start) {
			path.push_back(current);
		}
		visited[current] = true;
		degree = 0;
		edgeSum = 0;
		Graph::Node nextNode = lemon::INVALID;
		Graph::Node forbidden;
		if (path.empty()) {
			forbidden = exclude;
		} else if (path.size()==1) {
			forbidden = start;
		} else {
			forbidden = path[path.size()-2];
		}
		for (Graph::OutArcIt it(graph, current); it!=lemon::INVALID; ++it) {
			Graph::Node potential = graph.target(it);
			if (potential!=forbidden && (tolerance.negative(nextEdgeVal) || !tolerance.different(c[it], nextEdgeVal))) {
				nextNode = potential;
				if (tolerance.negative(nextEdgeVal)) {
					firstEdgeVal = c[it];
				}
				nextEdgeVal = c[it];
			}
			++degree;
			edgeSum += c[it];
		}
		nextEdgeVal = 1-nextEdgeVal;
		current = nextNode;
	} while (graph.valid(current) && !visited[current] && odd[current] && degree==2 &&
			 !tolerance.different(edgeSum, 1));
	if (degree==2 && graph.valid(current) && !tolerance.different(edgeSum, 1)) {
		path.push_back(current);
		visited[current] = true;
	}
	return path;
}