#include <utility>
#include <branch_and_cut.hpp>
#include <cmath>
#include <limits>
#include <ctime>

BranchAndCut::BranchAndCut(LinearProgram& program, const std::vector<CutGenerator*>& gens, size_t maxOpenSize) :
		problem(program), varCount(program.getVariableCount()), goal(program.getGoal()),
		currBest(static_cast<size_t>(varCount)), fractOpt(static_cast<size_t>(varCount)),
		objCoefficients(static_cast<size_t>(varCount)), maxOpenSize(maxOpenSize), generators(gens),
		constraintsAtStart(static_cast<size_t>(program.getConstraintCount())),
		defaultBounds(static_cast<size_t>(varCount)),
		currentBounds(static_cast<size_t>(varCount)),
		intTolerance(0.0001) {
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
long tolerantCeil(double in, lemon::Tolerance<double> tol) {
	long ceil = static_cast<long>(std::ceil(in));
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
long getNextWorse(double in, LinearProgram::Goal g, lemon::Tolerance<double> tol) {
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
void BranchAndCut::solveLP(LinearProgram::Solution& out) {
	CutGenerator::CutStatus solutionValid;
	double oldVal = 0;
	size_t slowIterations = 0;
	unsigned iterations = 0;
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
		if (!out.isValid() || !isBetter(getNextWorse(out.getValue(), goal, intTolerance), upperBound, goal)) {
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
		solutionValid = CutGenerator::valid;
		for (CutGenerator* gen:generators) {
			CutGenerator::CutStatus genStatus = gen->validate(problem, out.getVector());
			//a>b bedeutet, dass a eine Neuberechnung "dringender" macht als b
			if (genStatus>solutionValid) {
				solutionValid = genStatus;
			}
		}
		//sinceSlack0 vergrößern, da eventuell Constraints hinzugefügt wurden
		sinceSlack0.resize(problem.getConstraintCount()-constraintsAtStart, 0);
		if (slowIterations>6 && solutionValid==CutGenerator::maybe_recalc) {
			solutionValid = CutGenerator::valid;
		}
	} while (solutionValid!=CutGenerator::valid);
	const double maxRatio = 3;
	if (problem.getConstraintCount()>constraintsAtStart*maxRatio) {
		std::vector<int> toRemove(problem.getConstraintCount(), 0);
		size_t removeCount = 0;
		//Alle constraints, die seit 10 oder mehr Iterationen nicht mit Gleichheit erfüllt waren, werden entfernt
		for (size_t i = 0; i<sinceSlack0.size(); ++i) {
			if (sinceSlack0[i]>10) {
				toRemove[i+constraintsAtStart] = 1;
				++removeCount;
			}
		}
		if (!toRemove.empty()) {
			sinceSlack0.resize(sinceSlack0.size()-removeCount);
			std::fill(sinceSlack0.begin(), sinceSlack0.end(), 0);
			problem.removeSetConstraints(toRemove);
		}
	}
}

/**
 * @return eine optimale ganzzahlige Lösung des LP, die von alle Cut-Generatoren akzeptiert wird.
 */
std::vector<long> BranchAndCut::solve() {
	//Wurzel des Suchbaums zur offenen Menge hinzufügen
	BranchNode initNode{SystemBounds(this), 0, goal};
	open.insert(initNode);
	openSize += initNode.estimateSize();
	while (!open.empty()) {
		//besten offenen Knoten behandeln
		auto it = open.begin();
		BranchNode next = *it;
		open.erase(it);
		openSize -= next.estimateSize();
		//maxOpenSize==0->Es dürfen keine offenen Knoten gespeichert werden->DFS nutzen
		branchAndBound(next, maxOpenSize==0);
		std::cout << openSize << ", size: " << open.size() << std::endl;
		if (maxOpenSize>0) {
			/*
			 * Falls die offenen Menge "groß" wird, werden die Knoten mit dem höchsten Speicherverbrauch entfernt (mit
			 * DFS behandelt), bis die offenen Menge "klein" ist
			 */
			while (openSize>maxOpenSize*0.95) {
				auto maxIt = open.begin();
				size_t size = maxIt->estimateSize();
				for (it = open.begin(); it!=open.end(); ++it) {
					if (it->estimateSize()>size) {
						maxIt = it;
						size = it->estimateSize();
					}
				}
				next = *maxIt;
				std::cout << "Removing node of size " << next.estimateSize() << std::endl;
				openSize -= next.estimateSize();
				open.erase(maxIt);
				branchAndBound(next, true);
				std::cout << "Removed node, now at open size " << openSize << std::endl;
			}
		}
	}
	return currBest;
}

/**
 * Findet die beste ganzzahlige Lösung des LP's (mit Cut-Generatoren) unter den aktuellen Grenzen für die Variablen.
 * Falls diese Lösung besser als upperBound bzw.
 */
void BranchAndCut::branchAndBound(BranchNode& node, bool dfs) {
	setupBounds(node.bounds);
	solveLP(fractOpt);
	/*
	 * Ein Branch wird nur dann weiterentwickelt, wenn das LP lösbar ist und die Lösung besser ist als die obere Schranke.
	 * Da sowohl die Variablen als auch die Koeffizienten der Zielfunktion ganzzahlig sind, ist auch der Wert der
	 * Zielfunktion im Endergebnis ganzzahlig; sie hat also mindestens den Wert
	 * getNextWorse(fractOpt.getValue(), goal, intTolerance)
	 */
	if (fractOpt.isValid() && isBetter(getNextWorse(fractOpt.getValue(), goal, intTolerance), upperBound, goal)) {
		variable_id varToBound = -1;
		double optDist = 1;
		long varWeight = 0;
		const std::vector<double>& reducedCosts = fractOpt.getReducedCosts();
		size_t nonIntCount = 0;
		for (variable_id i = 0; i<varCount; ++i) {
			long rounded = std::lround(fractOpt[i]);
			double diff = rounded-fractOpt[i];
			if (intTolerance.nonZero(diff)) {
				/*
				 * Aktueller Wert der Variablen i ist nicht ganzzahlig
				 * Die Variable, die beschränkt wird, wird so gewählt, dass der Abstand zur nächsten ganzen Zahl maximal
				 * ist und (unter diesen Variablen) die Kosten der Variable betragsmäßig maximal sind.
				 */
				double distInt = std::abs(diff);
				//TODO Obj.-Koeff. oder reduzierte Kosten?
				long cost = std::abs(objCoefficients[i]);
				if (varToBound<0 || generalTolerance.less(optDist, distInt) ||
					(!generalTolerance.less(distInt, optDist) && cost>varWeight)) {
					optDist = distInt;
					varToBound = i;
					varWeight = cost;
				}
				++nonIntCount;
			} else if (diff==0) {
				LinearProgram::BoundType upper = LinearProgram::upper, lower = LinearProgram::lower;
				if (goal==LinearProgram::maximize) {
					std::swap(upper, lower);
				}
				/*
				 * Falls die Variable noch nicht komplett festgelegt ist, aber als Wert eine ihrer Schranken annimmt und
				 * die reduzierten Kosten bestimmte Bedingungen erfüllen, kann die Variable auf diesen Wert festgelegt
				 * werden
				 */
				if (node.bounds[i][upper]!=node.bounds[i][lower]) {
					double costDifference = upperBound-fractOpt.getValue();
					if (rounded==node.bounds[i][upper] && generalTolerance.less(reducedCosts[i], -costDifference)) {
						node.bounds.fix(i, upper);
					} else if (rounded==node.bounds[i][lower] &&
							   generalTolerance.less(costDifference, reducedCosts[i])) {
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
			std::vector<long> newOpt(static_cast<size_t>(varCount));
			long value = 0;
			for (variable_id i = 0; i<varCount; ++i) {
				newOpt[i] = std::lround(fractOpt[i]);
				value += newOpt[i]*objCoefficients[i];
			}
			if (value!=std::lround(fractOpt.getValue()))
				std::cerr << "Fractional: " << fractOpt.getValue() << ", integer: " << value << std::endl;
			setUpperBound(newOpt, value);
		} else {
			long rounded = std::lround(fractOpt[varToBound]);
			double diff = rounded-fractOpt[varToBound];
			long ceil;
			if (diff>0) {
				ceil = rounded;
			} else {
				ceil = rounded+1;
			}
			long floor = ceil-1;
			double fractVal = fractOpt.getValue();
			//Sofort bearbeiten, es sei denn, die offene Menge wird zu groß
			branch(varToBound, ceil, LinearProgram::lower, node.bounds, fractVal,
				   openSize<maxOpenSize, dfs);
			/*
			 * Später bearbeiten, es sei denn, die offene Menge ist klein und es gibt nur wenige nicht ganzzahlige
			 * Variablen
			 */
			branch(varToBound, floor, LinearProgram::upper, node.bounds, fractVal,
				   nonIntCount<10 && openSize<maxOpenSize, dfs);
		}
	}
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
void BranchAndCut::branch(variable_id variable, long val, LinearProgram::BoundType bound,
						  const SystemBounds& parent, double objValue, bool immediate, bool dfs) {
	BranchNode node{parent, objValue, goal};
	if (!node.bounds.bounds.count(variable)) {
		node.bounds.bounds[variable] = defaultBounds[variable];
	}
	node.bounds.bounds[variable][bound] = val;
	if (immediate || dfs) {
		branchAndBound(node, dfs);
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
bool BranchAndCut::isBetter(double a, double b, LinearProgram::Goal goal) {
	lemon::Tolerance<double> tolerance;
	if (goal==LinearProgram::maximize) {
		return tolerance.less(b, a);
	} else {
		return tolerance.less(a, b);
	}
}

/**
 * Setzt die obere Schranke für Branch and Bound
 * @param value Die Werte der Variablen
 * @param cost Die Kosten der angegebenen Lösung
 */
void BranchAndCut::setUpperBound(const std::vector<long>& value, long cost) {
	upperBound = cost;
	currBest = value;
	std::cout << "Setting upper bound as " << cost << std::endl;
	if (!open.empty()) {
		//TODO geht das schöner?
		auto it = open.end();
		--it;
		while (!open.empty() && !isBetter(getNextWorse(it->value, goal, intTolerance), cost, goal)) {
			openSize -= it->estimateSize();
			it = open.erase(it);
			if (!open.empty()) {
				--it;
			}
		}
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
			long newBound = boundsForVar[type];
			if (generalTolerance.different(currBound, newBound)) {
				problem.setBound(i, type, newBound);
				currentBounds[i][type] = newBound;
			}
		}
	}
}

long& BranchAndCut::VariableBounds::operator[](LinearProgram::BoundType b) {
	if (b==LinearProgram::lower) {
		return min;
	} else {
		return max;
	}
}

/**
 * Gibt eine (grobe) Schätzung für den Speicherverbrauch dieses Knotens zurück. Geht davon aus, dass std::vector<bool>
 * effizient implementiert wurde, d.h. dass in einem byte CHAR_BIT Elemente gespeichert werden
 */
size_t BranchAndCut::BranchNode::estimateSize() const {
	return sizeof(*this)+bounds.bounds.size()*sizeof(*bounds.bounds.begin())+
		   2*bounds.fixLower.size()/CHAR_BIT;
}

bool BranchAndCut::BranchNode::operator<(const BranchAndCut::BranchNode& other) const {
	return value<other.value;
}

BranchAndCut::SystemBounds::SystemBounds(BranchAndCut* owner)
		: owner(owner),
		  fixLower(static_cast<size_t>(owner->varCount), false),
		  fixUpper(static_cast<size_t>(owner->varCount), false) {}

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
	if (b==LinearProgram::lower) {
		fixLower[var] = true;
	} else {
		fixUpper[var] = true;
	}
}
