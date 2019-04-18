#include <utility>

#include <utility>
#include <branch_and_cut.hpp>
#include <cmath>
#include <limits>
#include <ctime>
#include <tsp_utils.hpp>

BranchAndCut::BranchAndCut(LinearProgram& program, std::vector<CutGenerator *>  gens, VariableRemover *remover,
						   size_t maxOpenSize) :
		problem(program), varCount(program.getVariableCount()), goal(program.getGoal()),
		currBest(static_cast<size_t>(varCount)), fractOpt(static_cast<size_t>(varCount), 0),
		objCoefficients(static_cast<size_t>(varCount)), maxOpenSize(maxOpenSize), generators(std::move(gens)),
		constraintsAtStart(static_cast<size_t>(program.getConstraintCount())),
		defaultBounds(static_cast<size_t>(varCount)), currentBounds(static_cast<size_t>(varCount)), intTolerance(0.01),
		remover(remover) {
	//obere Schranke auf schlechtesten möglichen Wert setzen
	if (goal==LinearProgram::minimize) {
		upperBound = std::numeric_limits<double>::max();
	} else {
		upperBound = -std::numeric_limits<double>::max();
	}
	//Koeffizienten in einem vector speichern, da direkter Zugriff über das LP langsam ist
	std::vector<double> obj = program.getObjective();
	for (size_t i = 0; i<obj.size(); ++i) {
		objCoefficients[i] = std::lround(obj[i]);
	}
	//Variablengrenzen auslesen und als Standardgrenzen speichern
	for (variable_id i = 0; i<varCount; ++i) {
		for (LinearProgram::BoundType type:{LinearProgram::lower, LinearProgram::upper}) {
			defaultBounds[i][type] = std::lround(problem.getBound(i, type));
		}
		currentBounds[i] = defaultBounds[i];
	}
}

/**
 * Berechnet die kleinste ganze Zahl, die mit der angegebenen Toleranz tol nicht kleiner als in ist
 */
value_t tolerantCeil(double in, lemon::Tolerance<double> tol) {
	auto ceil = static_cast<value_t>(std::ceil(in));
	if (tol.different(ceil-in, 1)) {
		return ceil;
	} else {
		return ceil-1;
	}
}

/**
 * Berechnet den besten möglichen ganzzahligen Wert, der schlechter als in ist. Besser/Schlechter bezieht sich auf g,
 * außerdem wird Ganzzahligkeit mit der Toleranz tol entschieden
 */
value_t getNextWorse(double in, LinearProgram::Goal g, lemon::Tolerance<double> tol) {
	if (g==LinearProgram::minimize) {
		return tolerantCeil(in, tol);
	} else {
		return -tolerantCeil(-in, tol);
	}
}

/**
 * Löst das LP wiederholt und fügt die von den Cut-Gereratoren erzeugten Schnitte hinzu, bis eine der folgenden
 * Bedingungen erfüllt ist:<br>
 * 1. Das LP ist unzulässig.<br>
 * 2. Die Lösung wird von allen Cut-Generatoren als gültig akzeptiert.<br>
 * 3. Die Lösung ist schlechter als die aktuelle obere Schranke.<br>
 * 4. Der Wert der Lösung hat sich in den letzten Iterationen nicht stark geändert und es wurde kein Schnitt hinzu-
 * gefügt, der eine Neuberechnung erzwingt.<br>
 * Am Ende werden Constraints entfernt, die seit mindestens 10 Iterationen nicht mehr mit Gleichheit erfüllt waren.
 * @param out Ein Solution-Objekt der korrekten Größe. Dient als Ausgabe.
 */
void BranchAndCut::solveLP(LinearProgram::Solution& out, bool isRoot) {
	CutGenerator::CutStatus solutionStatus;
	double oldVal = 0;
	const size_t maxSlow = 6;
	size_t slowIterations = 0;
	unsigned iterations = 0;
	double lastCleanup = 0;
	do {
		++iterations;
		if (iterations%64==0) {
			std::cout << "Currently at " << iterations << " iterations, using " << problem.getConstraintCount()
					  << " constraints!" << std::endl;
		}
		problem.solve(out);
		countSolutionSlack(out);
		/*
		 * Abbrechen, falls das LP unzulässig ist oder die optimale fraktionale Lösung nicht mehr besser als die beste
		 * bekannte ganzzahlige ist
		 */
		if (!out.isValid() || !isBetter(getNextWorse(out.getValue(), goal, intTolerance), upperBound)) {
			break;
		}
		/*
		 * Als langsame Iteration zählen solche, in denen die relative Änderung des Wertes der Lösung weniger
		 * als 2.5*10^-5 ist, siehe "A branch-and-cut algorithm for the resolution of large-scale symmetric traveling
		 * salesman problems", Padberg und Rinaldi 1991
		 */
		if (std::abs(oldVal-out.getValue())<oldVal*2.5e-5) {
			++slowIterations;
		} else {
			slowIterations = 0;
		}
		oldVal = std::abs(out.getValue());
		solutionStatus = readdRemovedConstraints(out);
		if (solutionStatus==CutGenerator::valid) {
			for (CutGenerator* gen:generators) {
				//Kurz
				if (slowIterations > maxSlow - 1 && solutionStatus == CutGenerator::maybe_recalc) {
					solutionStatus = CutGenerator::valid;
				}
				CutGenerator::CutStatus genStatus = gen->validate(problem, out.getVector(), solutionStatus);
				//a>b bedeutet, dass a eine Neuberechnung "dringender" macht als b
				if (genStatus>solutionStatus) {
					solutionStatus = genStatus;
				}
			}
		}
		//sinceSlack0 vergrößern, da Constraints hinzugefügt wurden
		if (solutionStatus!=CutGenerator::valid) {
			sinceSlack0.resize(problem.getConstraintCount()-constraintsAtStart, 0);
		}
		if (slowIterations > maxSlow && solutionStatus == CutGenerator::maybe_recalc) {
			solutionStatus = CutGenerator::valid;
		}
		if (iterations == 1) {
			lastCleanup = out.getValue();
		}
		if (iterations==1 || (iterations%8==0 && isBetter(lastCleanup, out.getValue()))) {
			/*
			 * Entfernen "alter" Constraints maximal alle 8 Iterationen und nur, wenn der Wert der aktuellen Lösung
			 * "schlechter" (also eine bessere untere Schranke) als beim letzten erfolgreichen Entfernen ist
			 */
			if (cleanupOldConstraints()) {
				lastCleanup = out.getValue();
			}
			if (remover!=nullptr && isRoot) {
				std::vector<variable_id> toRemove = remover->removeOnRootSolution(out);
				removeVariables(toRemove);
			}
		}
	} while (solutionStatus!=CutGenerator::valid);
}

/**
 * @return eine optimale ganzzahlige Lösung des LP, die von alle Cut-Generatoren akzeptiert wird.
 */
std::vector<value_t> BranchAndCut::solve() {
	handledNodes = 0;
	//Wurzel des Suchbaums zur offenen Menge hinzufügen
	BranchNode initNode{SystemBounds(this), 0, goal};
	open.insert(initNode);
	openSize += initNode.estimateSize();
	size_t iterations = 0;
	bool isRoot = true;
	while (!open.empty()) {
		//besten offenen Knoten behandeln
		auto it = open.begin();
		BranchNode next = *it;
		open.erase(it);
		openSize -= next.estimateSize();
		//maxOpenSize==0->Es dürfen keine offenen Knoten gespeichert werden->DFS nutzen
		branchAndBound(next, maxOpenSize == 0, isRoot);
		isRoot = false;
		++iterations;
		if (iterations%16==0) {
			std::cout << "At iteration " << iterations << ", cost is " << next.value << std::endl;
		}
		if (maxOpenSize>0) {
			/*
			 * Falls die offenen Menge "groß" wird, werden die Knoten entfernt (mit DFS behandelt), bis die offenen Menge
			 * "klein" ist. Dabei werden Knoten gewählt, die möglichst viel Speicher verbrauchen, unter diesen solche mit
			 * möglichst vielen festgelegten Variablen und unter diesen solche mit der kleinsten unteren Schranke. Die
			 * erste Bedingung führt dazu, dass der Speicherverbrauch möglicht schnell klein wird; die zweite Bedingung
			 * soll bewirken, dass die Knoten oft sehr schnell bearbeitet werden können; die letzte Bedingung soll dazu
			 * führen, dass möglichst wenige Knoten bearbeitet werden, die durch eine bessere obere Schranke
			 * ausgeschlossen werden könnten.
			 */
			while (openSize>maxOpenSize*0.95) {
				std::cout << "Open set has size " << openSize << ", looking for node to remove" << std::endl;
				auto maxIt = open.begin();
				size_t maxSize = maxIt->estimateSize();
				size_t maxFixed = 0;
				for (it = open.begin(); it!=open.end(); ++it) {
					size_t estSize = it->estimateSize();
					size_t fixedCount = it->getFixedCount();
					if (estSize>maxSize || (estSize==maxSize && fixedCount>maxFixed)) {
						maxIt = it;
						maxSize = estSize;
						maxFixed = fixedCount;
					}
				}
				next = *maxIt;
				std::cout << "Removing node of size " << next.estimateSize() << std::endl;
				openSize -= next.estimateSize();
				open.erase(maxIt);
				branchAndBound(next, true, false);
				std::cout << "Removed node, now at open size " << openSize << std::endl;
			}
		}
	}
	std::cout << "Found solution in " << handledNodes << " nodes" << std::endl;
	return currBest;
}

/**
 * Findet die beste ganzzahlige Lösung des LP's (mit Cut-Generatoren) unter den aktuellen Grenzen für die Variablen.
 * Falls diese Lösung besser als upperBound bzw.
 */
void BranchAndCut::branchAndBound(BranchNode& node, bool dfs, bool isRoot) {
	BranchNode* prevNode = currentNode;
	currentNode = &node;
	++handledNodes;
	setupBounds(node.bounds);
	solveLP(fractOpt, isRoot);
	/*
	 * Ein Branch wird nur dann weiterentwickelt, wenn das LP lösbar ist und die Lösung besser ist als die obere Schranke.
	 * Da sowohl die Variablen als auch die Koeffizienten der Zielfunktion ganzzahlig sind, ist auch der Wert der
	 * Zielfunktion im Endergebnis ganzzahlig; sie hat also mindestens den Wert
	 * getNextWorse(fractOpt.getValue(), goal, intTolerance)
	 */
	if (fractOpt.isValid() && isBetter(getNextWorse(fractOpt.getValue(), goal, intTolerance), upperBound)) {
		variable_id varToBound = -1;
		double optDist = 1;
		coeff_t varWeight = 0;
		const std::vector<double>& reducedCosts = fractOpt.getReducedCosts();
		size_t nonIntCount = 0;
		for (variable_id i = 0; i<varCount; ++i) {
			value_t rounded = std::lround(fractOpt[i]);
			double diff = rounded-fractOpt[i];
			if (intTolerance.nonZero(diff)) {
				/*
				 * Aktueller Wert der Variablen i ist nicht ganzzahlig
				 * Die Variable, die beschränkt wird, wird so gewählt, dass der Abstand zur nächsten ganzen Zahl maximal
				 * ist und (unter diesen Variablen) die Kosten der Variable betragsmäßig maximal sind.
				 */
				double distInt = std::abs(diff);
				coeff_t cost = std::abs(objCoefficients[i]);
				if (varToBound<0 || generalTolerance.less(optDist, distInt) ||
					(!generalTolerance.less(distInt, optDist) && cost>varWeight)) {
					optDist = distInt;
					varToBound = i;
					varWeight = cost;
				}
				++nonIntCount;
			} else {
				LinearProgram::BoundType upper = LinearProgram::upper, lower = LinearProgram::lower;
				//TODO stimmt das so für max.?
				if (goal==LinearProgram::maximize) {
					std::swap(upper, lower);
				}
				/*
				 * Falls die Variable noch nicht komplett festgelegt ist, aber als Wert eine ihrer Schranken annimmt und
				 * die reduzierten Kosten bestimmte Bedingungen erfüllen, kann die Variable auf diesen Wert festgelegt
				 * werden
				 */
				if (node.bounds[i][upper]!=node.bounds[i][lower]) {
					//TODO upperBound ist für das ganze IP, nicht für diesen Zweig. Geht das trotzdem?
					//Nach Intuition ja, für ein richtiges Argument brauche ich mehr LGO. -1 ist auch nach Intuition
					double minDifference = upperBound-1-fractOpt.getValue();
					if (rounded==node.bounds[i][upper] && intTolerance.less(reducedCosts[i], -minDifference)) {
						node.bounds.fix(i, upper);
					} else if (rounded==node.bounds[i][lower] && intTolerance.less(minDifference, reducedCosts[i])) {
						node.bounds.fix(i, lower);
					}
				}
			}
		}
		if (varToBound<0) {
			/*
			 * Es wurde keine nicht-ganzzahlige Variable gefunden. Außerdem ist der Wert der Lösung besser als die
			 * bisherige obere Schranke; die Lösung wird also als neue obere Schranke gesetzt.
			 */
			std::vector<value_t> newOpt(static_cast<size_t>(varCount));
			coeff_t value = 0;
			for (variable_id i = 0; i<varCount; ++i) {
				newOpt[i] = std::lround(fractOpt[i]);
				value += newOpt[i]*objCoefficients[i];
			}
			if (value!=std::lround(fractOpt.getValue()))
				std::cerr << "Fractional: " << fractOpt.getValue() << ", integer: " << value << std::endl;
			setUpperBound(newOpt, value);
		} else {
			value_t rounded = std::lround(fractOpt[varToBound]);
			double diff = rounded-fractOpt[varToBound];
			value_t ceil;
			if (diff>0) {
				ceil = rounded;
			} else {
				ceil = rounded+1;
			}
			value_t floor = ceil-1;
			double fractVal = fractOpt.getValue();
			/*
			 * Später bearbeiten, es sei denn, die offene Menge ist klein und es gibt nur wenige nicht ganzzahlige
			 * Variablen
			 */
			branch(varToBound, floor, LinearProgram::upper, node.bounds, fractVal,
				   nonIntCount < 10 && openSize < maxOpenSize && !isRoot, dfs);
			/*
			 * Sofort bearbeiten, es sei denn, die offene Menge wird zu groß. So wird die offene Menge relativ klein
			 * gehalten, es werden schneller obere Schranken gefunden und die LPs können schneller gelöst werden, da
			 * das duale Simplexverfahren genutzt wird und das LP nach dem hinzufügen neuer Constraints schnell neu
			 * gelöst werden kann (TODO genauer beschreiben, wenn das in LGO besprochen wurde)
			 * TODO warum wird das beim vertauschen von lower und upper deutlich langsamer?
			 */
			branch(varToBound, ceil, LinearProgram::lower, node.bounds, fractVal,
				   openSize<maxOpenSize && !isRoot, dfs);
		}
	}
	currentNode = prevNode;
}

/**
 * Beschränkt den zulässgen Bereich einer Variablen wie durch die Parameter beschrieben, ruft dann branchAndBound auf
 * und setzt den zulässigen Bereich wieder auf den ursprünglichen Wert.
 * @param variable Die zu beschränkende Variable
 * @param val Der neue Wert der Beschränkung
 * @param bound Die "Richtung", in der die Variable beschränkt werden soll
 * @param parent Die sonstigen Variablenbeschränkungen
 * @param objValue Der optimale (fraktionale) Wert der Zielfunktion im Vorgängerknoten (d.h. mit parent als Beschränkung)
 * @param immediate true, falls der Knoten sofort behandelt werden soll. False, falls der Knoten zur Menge der offenen
 * Knoten hinzugefügt werden soll
 * @param dfs Falls true, werden alle Kinder des entstehenden Knoten sofort behandelt (immediate = true)
 */
void BranchAndCut::branch(variable_id variable, value_t val, LinearProgram::BoundType bound,
						  const SystemBounds& parent, double objValue, bool immediate, bool dfs) {
	BranchNode node{parent, objValue, goal};
	node.bounds.setBound(variable, bound, val);
	if (immediate || dfs) {
		branchAndBound(node, dfs, false);
	} else {
		open.insert(node);
		openSize += node.estimateSize();
	}
}

/**
 * Erhöht die Zähler für die Constraints, die nicht mit Gleichheit erfüllt sind und setzt die Zähler für die anderen
 * zurück.
 */
void BranchAndCut::countSolutionSlack(const LinearProgram::Solution& sol) {
	const std::vector<double>& slack = sol.getSlack();
	for (size_t constraint = constraintsAtStart; constraint<slack.size(); ++constraint) {
		if (generalTolerance.nonZero(slack[constraint])) {
			++sinceSlack0[constraint-constraintsAtStart];
		} else {
			sinceSlack0[constraint-constraintsAtStart] = 0;
		}
	}
}

/**
 * @param a ein Wert der Zielfunktion
 * @param b ein Wert der Zielfunktion
 * @return true, falls a ein streng besserer Wert der Zielfunktion als b ist
 */
bool BranchAndCut::isBetter(double a, double b) {
	if (goal==LinearProgram::maximize) {
		return generalTolerance.less(b, a);
	} else {
		return generalTolerance.less(a, b);
	}
}

/**
 * Setzt die obere Schranke für Branch and Bound
 * @param value Die Werte der Variablen
 * @param cost Die Kosten der angegebenen Lösung
 */
void BranchAndCut::setUpperBound(const std::vector<value_t>& value, coeff_t cost) {
	upperBound = cost;
	if (value.size() != varCount) {
		std::cout << "Setting upper bound as " << upperBound << " without an example vector" << std::endl;
	} else {
		currBest = value;
		std::cout << "Setting upper bound as " << upperBound << std::endl;
	}
	if (!open.empty()) {
		size_t removed = 1;
		auto firstRemove = open.end();
		--firstRemove;
		while (firstRemove!=open.begin() && !isBetter(getNextWorse(firstRemove->value, goal, intTolerance), cost)) {
			openSize -= firstRemove->estimateSize();
			--firstRemove;
			++removed;
		}
		if (isBetter(getNextWorse(firstRemove->value, goal, intTolerance), cost)) {
			++firstRemove;
			--removed;
		}
		open.erase(firstRemove, open.end());
		std::cout << "Removed " << removed << " nodes from the open set" << std::endl;
	}
	if (value.size() == varCount && remover != nullptr) {
		std::vector<variable_id> toRemove = remover->removeOnUpperBound(currBest);
		removeVariables(toRemove);
	}
}

/**
 * Setzt die Variablenbeschränkungen des LPs auf die angegebenen Werte
 */
void BranchAndCut::setupBounds(const SystemBounds& bounds) {
	for (variable_id i = 0; i<varCount; ++i) {
		VariableBounds boundsForVar = bounds[i];
		for (LinearProgram::BoundType type:{LinearProgram::lower, LinearProgram::upper}) {
			double currBound = currentBounds[i][type];
			value_t newBound = boundsForVar[type];
			if (generalTolerance.different(currBound, newBound)) {
				problem.setBound(i, type, newBound);
				currentBounds[i][type] = newBound;
			}
		}
	}
}

/**
 * Fügt vor kurzem entfernte Constraints zum LP hinzu, die von der übergebenen Lösung verletzt werden
 * @return valid, falls keine Constraints hinzugefügt wurden, sonst recalc
 */
CutGenerator::CutStatus BranchAndCut::readdRemovedConstraints(const LinearProgram::Solution& sol) {
	CutGenerator::CutStatus ret = CutGenerator::valid;
	size_t i = 0;
	while (i<recentlyRemoved.size()) {
		const LinearProgram::Constraint& constr = recentlyRemoved[i];
		double lhs = constr.evalLHS(sol.getVector());
		//TODO intTolerance, general oder was anderes?
		if (!constr.isValidLHS(lhs, intTolerance)) {
			ret = CutGenerator::recalc;
			problem.addConstraint(constr);
			recentlyRemoved[i] = std::move(recentlyRemoved.back());
			recentlyRemoved.pop_back();
		} else {
			++i;
		}
	}
	return ret;
}

bool BranchAndCut::cleanupOldConstraints() {
	const double maxRatio = 2;
	const size_t minNonEqual = 10;
	const size_t minRemove = constraintsAtStart / 10;
	if (problem.getConstraintCount() > constraintsAtStart * maxRatio) {
		std::vector<int> toRemove(problem.getConstraintCount(), 0);
		size_t removeCount = 0;
		std::vector<LinearProgram::Constraint> removed;
		//Alle constraints, die seit 10 oder mehr Iterationen nicht mit Gleichheit erfüllt waren, werden entfernt
		for (size_t i = 0; i < sinceSlack0.size(); ++i) {
			if (sinceSlack0[i] > minNonEqual) {
				size_t indexInLP = i + constraintsAtStart;
				toRemove[indexInLP] = 1;
				++removeCount;
				removed.push_back(problem.getConstraint(static_cast<int>(indexInLP)));
			}
		}
		if (removeCount > minRemove) {
			sinceSlack0.resize(sinceSlack0.size() - removeCount);
			std::fill(sinceSlack0.begin(), sinceSlack0.end(), 0);
			problem.removeSetConstraints(toRemove);
			std::cout << "Removed " << removeCount << " old constraints" << std::endl;
			recentlyRemoved = removed;
			return true;
		}
	}
	return false;

}

void BranchAndCut::removeVariables(const std::vector<variable_id>& toRemove) {
	if (!toRemove.empty()) {
		std::vector<int> variableMap(varCount);
		for (variable_id r:toRemove) {
			variableMap[r] = 1;
		}
		problem.removeSetVariables(variableMap);
		variable_id newVarCount = varCount - toRemove.size();
		std::multiset<BranchNode> newOpen;
		for (const BranchNode& old:open) {
			SystemBounds newBounds(old.bounds, variableMap, newVarCount);
			BranchNode newNode{newBounds, old.value, old.goal};
			bool add = true;
			for (const BranchNode& newEle:newOpen) {
				if (newNode == newEle) {
					add = false;
					break;
				}
			}
			if (add) {
				newOpen.insert(newNode);
			}
		}
		if (currentNode!=nullptr) {
			SystemBounds newBounds(currentNode->bounds, variableMap, newVarCount);
			*currentNode = {newBounds, currentNode->value, currentNode->goal};
		}
		tsp_util::eraseEntries(defaultBounds, toRemove);
		tsp_util::eraseEntries(currentBounds, toRemove);
		tsp_util::eraseEntries(objCoefficients, toRemove);
		for (LinearProgram::Constraint& c:recentlyRemoved) {
			c.deleteVariables(variableMap);
		}
		varCount = newVarCount;
		open = newOpen;
		fractOpt.removeVariables(toRemove);
		assert(varCount == defaultBounds.size() && varCount == currentBounds.size() &&
			   varCount == objCoefficients.size());
		std::cout << "Removed " << toRemove.size() << " variables, " << varCount << " remaining" << std::endl;
	}
}

value_t& BranchAndCut::VariableBounds::operator[](LinearProgram::BoundType b) {
	if (b==LinearProgram::lower) {
		return min;
	} else {
		return max;
	}
}

bool BranchAndCut::VariableBounds::isFixed() const {
	return max==min;
}

bool BranchAndCut::VariableBounds::operator==(const BranchAndCut::VariableBounds& other) const {
	return min==other.min && max==other.max;
}

/**
 * Gibt eine (grobe) Schätzung für den Speicherverbrauch dieses Knotens zurück. Geht davon aus, dass std::vector<bool>
 * effizient implementiert wurde, d.h. dass in einem byte CHAR_BIT Elemente gespeichert werden
 */
size_t BranchAndCut::BranchNode::estimateSize() const {
	return sizeof(*this)+bounds.estimateSize()-sizeof(bounds);
}

size_t BranchAndCut::BranchNode::getFixedCount() const {
	return bounds.getFixedCount();
}

bool BranchAndCut::BranchNode::operator<(const BranchAndCut::BranchNode& other) const {
	return value<other.value;
}

bool BranchAndCut::BranchNode::operator==(const BranchAndCut::BranchNode& other) const {
	return bounds==other.bounds;
}

BranchAndCut::SystemBounds::SystemBounds(BranchAndCut* owner)
		: owner(owner),
		  fixLower(static_cast<size_t>(owner->varCount), false),
		  fixUpper(static_cast<size_t>(owner->varCount), false),
		  fixedCount(0) {}

BranchAndCut::SystemBounds::SystemBounds(const BranchAndCut::SystemBounds& old, const std::vector<variable_id>& idMap,
										 variable_id newVarCount) :
		owner(old.owner), fixedCount(0) {
	fixLower.resize(newVarCount);
	fixUpper.resize(newVarCount);
	for (size_t oldId = 0; oldId < old.fixLower.size(); ++oldId) {
		variable_id newId = idMap[oldId];
		if (newId >= 0) {
			fixUpper[newId] = old.fixUpper[oldId];
			fixLower[newId] = old.fixLower[oldId];
			if (old.bounds.count(oldId) > 0) {
				bounds[newId] = old.bounds.at(oldId);
			}
		}
	}
}

BranchAndCut::VariableBounds BranchAndCut::SystemBounds::operator[](variable_id id) const {
	VariableBounds basic{};
	if (bounds.count(id)) {
		basic = bounds.at(id);
	} else {
		basic = owner->defaultBounds[id];
	}
	if (fixUpper[id]) {
		return {basic.max, basic.max};
	} else if (fixLower[id]) {
		return {basic.min, basic.min};
	} else {
		return basic;
	}
}

/**
 * Setzt den Wert der Variablen var auf den durch b angegebenen "Extremwert" des aktuellen Bereichs fest
 */
void BranchAndCut::SystemBounds::fix(variable_id var, LinearProgram::BoundType b) {
	assert(!(*this)[var].isFixed());
	if (b==LinearProgram::lower) {
		fixLower[var] = true;
	} else {
		fixUpper[var] = true;
	}
	++fixedCount;
}

void BranchAndCut::SystemBounds::setBound(variable_id var, LinearProgram::BoundType b, value_t newBound) {
	VariableBounds varBounds = (*this)[var];
	bool wasFixed = varBounds.isFixed();
	varBounds[b] = newBound;
	if (varBounds.isFixed()) {
		if (!wasFixed) {
			fix(var, b==LinearProgram::lower ? LinearProgram::upper : LinearProgram::lower);
		}
	} else {
		bounds[var] = varBounds;
		if (wasFixed) {
			--fixedCount;
		}
	}
}

bool BranchAndCut::SystemBounds::operator==(const BranchAndCut::SystemBounds& other) const {
	for (variable_id i = 0; i<fixLower.size(); ++i) {
		if (!((*this)[i]==other[i])) {
			return false;
		}
	}
	return true;
}

size_t BranchAndCut::SystemBounds::getFixedCount() const {
	return fixedCount;
}

size_t BranchAndCut::SystemBounds::estimateSize() const {
	return bounds.size()*sizeof(*bounds.begin())+
		   2*fixLower.size()/CHAR_BIT+sizeof(*this);
}
