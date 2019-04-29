#include <utility>
#include <branch_and_cut.hpp>
#include <cmath>
#include <limits>
#include <tsp_utils.hpp>
#include <numeric>
#include <cassert>
#include <lemon/tolerance.h>
#include <algorithm>
#include <iostream>
#include <set>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdlib>
#include <cut_generator.hpp>
#include <linear_program.hpp>
#include <relative_tolerance.hpp>
#include <variable_remover.hpp>

using std::size_t;

BranchAndCut::BranchAndCut(LinearProgram& program, std::vector<CutGenerator *> gens, VariableRemover *remover,
						   bool dfs) :
		problem(program), varCount(program.getVariableCount()),
		currBest(static_cast<size_t>(varCount)), fractOpt(static_cast<size_t>(varCount), 0),
		objCoefficients(static_cast<size_t>(varCount)), currentBounds(static_cast<size_t>(varCount)),
		dfs(dfs), generators(std::move(gens)), remover(remover),
		constraintsAtStart(static_cast<size_t>(program.getConstraintCount())),
		intTolerance(0.01), correspondingOrigVar(varCount) {
	//Prüfen, dass das LP die Bedingungen für diese Klasse erfüllt
	if (program.getGoal() != LinearProgram::minimize) {
		throw std::runtime_error("Problem is not a minimization problem!");
	}
	for (variable_id i = 0; i < varCount; ++i) {
		double lower = problem.getBound(i, LinearProgram::lower);
		if (lower != 0) {
			throw std::runtime_error("Variable " + std::to_string(i) + " has lower bound " + std::to_string(lower)
									 + ", expected 0");
		}
		double upper = problem.getBound(i, LinearProgram::upper);
		if (upper != 1) {
			throw std::runtime_error("Variable " + std::to_string(i) + " has upper bound " + std::to_string(upper)
									 + ", expected 1");
		}
	}
	//obere Schranke auf schlechtesten möglichen Wert setzen
	upperBound = std::numeric_limits<double>::max();
	//Koeffizienten in einem vector speichern, da direkter Zugriff über das LP langsam ist
	std::vector<double> obj = program.getObjective();
	for (size_t i = 0; i < obj.size(); ++i) {
		double realCoeff = obj[i];
		objCoefficients[i] = std::lround(realCoeff);
		if (realCoeff != objCoefficients[i]) {
			throw std::runtime_error("Objective coefficient for " + std::to_string(i) + " is noninteger: " +
									 "Difference to closest integer is " +
									 std::to_string(objCoefficients[i] - realCoeff));
		}
	}
	//Am Anfang sind ursprüngliche und aktuelle ID's gleich
	std::iota(correspondingOrigVar.begin(), correspondingOrigVar.end(), 0);
}

/**
 * Berechnet die kleinste ganze Zahl, die mit der angegebenen Toleranz tol nicht kleiner als in ist
 */
obj_t tolerantCeil(double in, lemon::Tolerance<double> tol) {
	auto ceil = static_cast<obj_t>(std::ceil(in));
	if (tol.different(ceil - in, 1)) {
		return ceil;
	} else {
		return ceil - 1;
	}
}

/**
 * Berechnet den besten möglichen ganzzahligen Wert, der schlechter als in ist. Besser/Schlechter bezieht sich auf g,
 * außerdem wird Ganzzahligkeit mit der Toleranz tol entschieden
 */
obj_t getNextWorse(double in, lemon::Tolerance<double> tol) {
	return tolerantCeil(in, tol);
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
		if (iterations % 64 == 0) {
			std::cout << "Currently at " << iterations << " cutting iterations, using " << problem.getConstraintCount()
					  << " constraints!" << std::endl;
		}
		problem.solve(out);
		countSolutionSlack(out);
		/*
		 * Abbrechen, falls das LP unzulässig ist oder die optimale fraktionale Lösung nicht mehr besser als die beste
		 * bekannte ganzzahlige ist
		 */
		if (!out.isValid() || !isBetter(getNextWorse(out.getValue(), intTolerance), upperBound)) {
			break;
		}
		/*
		 * Als langsame Iteration zählen solche, in denen die relative Änderung des Wertes der Lösung weniger
		 * als 2.5*10^-5 ist, siehe "A branch-and-cut algorithm for the resolution of large-scale symmetric traveling
		 * salesman problems", Padberg und Rinaldi 1991
		 */
		if (std::abs(oldVal - out.getValue()) < oldVal * 2.5e-5) {
			++slowIterations;
		} else {
			slowIterations = 0;
		}
		oldVal = std::abs(out.getValue());
		solutionStatus = readdRemovedConstraints(out);
		if (solutionStatus == CutGenerator::valid) {
			for (CutGenerator *gen:generators) {
				//Kurz
				if (slowIterations > maxSlow - 1 && solutionStatus == CutGenerator::maybe_recalc) {
					solutionStatus = CutGenerator::valid;
				}
				CutGenerator::CutStatus genStatus = gen->validate(problem, out.getVector(), solutionStatus);
				//a>b bedeutet, dass a eine Neuberechnung "dringender" macht als b
				if (genStatus > solutionStatus) {
					solutionStatus = genStatus;
				}
			}
		}
		//sinceSlack0 vergrößern, da Constraints hinzugefügt wurden
		if (solutionStatus != CutGenerator::valid) {
			sinceSlack0.resize(problem.getConstraintCount() - constraintsAtStart, 0);
		}
		if (slowIterations > maxSlow && solutionStatus == CutGenerator::maybe_recalc) {
			solutionStatus = CutGenerator::valid;
		}
		if (iterations == 1) {
			lastCleanup = out.getValue();
		}
		if (iterations == 1 || (iterations % 8 == 0 && isBetter(lastCleanup, out.getValue()))) {
			/*
			 * Entfernen "alter" Constraints maximal alle 8 Iterationen und nur, wenn der Wert der aktuellen Lösung
			 * "schlechter" (also eine bessere untere Schranke) als beim letzten erfolgreichen Entfernen ist
			 */
			if (cleanupOldConstraints()) {
				lastCleanup = out.getValue();
			}
			if (remover != nullptr && isRoot) {
				std::vector<variable_id> toRemove = remover->removeOnRootSolution(out);
				removeVariables(toRemove);
			}
		}
	} while (solutionStatus != CutGenerator::valid);
}

/**
 * @return eine optimale ganzzahlige Lösung des LP, die von alle Cut-Generatoren akzeptiert wird.
 */
std::vector<value_t> BranchAndCut::solve() {
	handledNodes = 0;
	//Wurzel des Suchbaums zur offenen Menge hinzufügen
	BranchNode initNode{SystemBounds(varCount), 0};
	open.insert(initNode);
	size_t iterations = 0;
	bool isRoot = true;
	while (!open.empty()) {
		//besten offenen Knoten behandeln
		auto it = open.begin();
		BranchNode next = *it;
		open.erase(it);
		branchAndBound(next, isRoot);
		isRoot = false;
		++iterations;
		if (iterations % 16 == 0) {
			std::cout << "At iteration " << iterations << ", cost is " << next.value << ", open set contains " <<
					  open.size() << " nodes" << std::endl;
		}
	}
	std::cout << "Found solution in " << handledNodes << " nodes" << std::endl;
	return currBest;
}

/**
 * Findet die beste ganzzahlige Lösung des LP's (mit Cut-Generatoren) unter den aktuellen Grenzen für die Variablen.
 * Falls diese Lösung besser als upperBound bzw.
 */
void BranchAndCut::branchAndBound(BranchNode& node, bool isRoot) {
	currStack.push_back(&node);
	++handledNodes;
	if (handledNodes % 16 == 0) {
		std::cout << "Handling search node " << handledNodes << std::endl;
	}
	setupBounds(node.bounds);
	solveLP(fractOpt, isRoot);
	/*
	 * Ein Branch wird nur dann weiterentwickelt, wenn das LP lösbar ist und die Lösung besser ist als die obere Schranke.
	 * Da sowohl die Variablen als auch die Koeffizienten der Zielfunktion ganzzahlig sind, ist auch der Wert der
	 * Zielfunktion im Endergebnis ganzzahlig; sie hat also mindestens den Wert
	 * getNextWorse(fractOpt.getValue(), goal, intTolerance)
	 */
	if (fractOpt.isValid() && isBetter(getNextWorse(fractOpt.getValue(), intTolerance), upperBound)) {
		//Variable, die zum Branchen genutzt werden soll
		variable_id varToBound = -1;
		//Abstand des Wertes zu 0.5
		double optDist = 1;
		//Koeffizient der Variablen in der Zielfunktion
		coeff_t varWeight = 0;
		const std::vector<double>& reducedCosts = fractOpt.getReducedCosts();
		size_t nonIntCount = 0;
		for (variable_id i = 0; i < varCount; ++i) {
			value_t rounded = std::lround(fractOpt[i]);
			double diff = rounded - fractOpt[i];
			if (intTolerance.nonZero(diff)) {
				/*
				 * Aktueller Wert der Variablen i ist nicht ganzzahlig
				 * Die Variable, die beschränkt wird, wird so gewählt, dass der Abstand zur nächsten ganzen Zahl maximal
				 * ist und (unter diesen Variablen) die Kosten der Variable betragsmäßig maximal sind.
				 */
				double distToInt = std::abs(diff);
				coeff_t cost = std::abs(objCoefficients[i]);
				if (varToBound < 0 || generalTolerance.less(optDist, distToInt) ||
					(!generalTolerance.less(distToInt, optDist) && cost > varWeight)) {
					optDist = distToInt;
					varToBound = i;
					varWeight = cost;
				}
				++nonIntCount;
			} else {
				/*
				 * Falls die Variable noch nicht komplett festgelegt ist, aber als Wert eine ihrer Schranken annimmt und
				 * die reduzierten Kosten bestimmte Bedingungen erfüllen, kann die Variable auf diesen Wert festgelegt
				 * werden
				 */
				if (!node.bounds[i].isFixed()) {
					double minDifference = upperBound - 1 - fractOpt.getValue();
					if (rounded == 1 && intTolerance.less(reducedCosts[i], -minDifference)) {
						node.bounds.fix(i, LinearProgram::upper);
					} else if (rounded == 0 && intTolerance.less(minDifference, reducedCosts[i])) {
						node.bounds.fix(i, LinearProgram::lower);
					}
				}
			}
		}
		if (varToBound < 0) {
			/*
			 * Es wurde keine nicht-ganzzahlige Variable gefunden. Außerdem ist der Wert der Lösung besser als die
			 * bisherige obere Schranke; die Lösung wird also als neue obere Schranke gesetzt.
			 */
			std::vector<value_t> newOpt(static_cast<size_t>(varCount));
			obj_t value = 0;
			for (variable_id i = 0; i < varCount; ++i) {
				newOpt[i] = std::lround(fractOpt[i]);
				value += newOpt[i] * objCoefficients[i];
			}
			if (value != std::lround(fractOpt.getValue())) {
				std::cerr << "Fractional: " << fractOpt.getValue() << ", integer: " << value << std::endl;
			}
			setUpperBound(newOpt, value);
		} else {
			double fractVal = fractOpt.getValue();
			variable_id origVar = correspondingOrigVar[varToBound];
			/*
			 * Später bearbeiten, es sei denn, es gibt nur wenige nicht ganzzahlige Variablen (dann scheint es
			 * wahrscheinlich, dass hier schnell eine neue obere Schranke gefunden werden kann)
			 */
			branch(varToBound, 0, LinearProgram::upper, node.bounds, fractVal, nonIntCount < 10);
			//Wurden Variablen entfernt, so dass sich die ID der aktuelleb Variable geändert hat?
			if (varToBound >= varCount || correspondingOrigVar[varToBound] != origVar) {
				//Neue ID finden oder feststellen, dass die Variable entfernt wurde (invalid_variable)
				varToBound = LinearProgram::invalid_variable;
				for (variable_id i = 0; i < varCount; ++i) {
					if (correspondingOrigVar[i] == origVar) {
						varToBound = i;
						break;
					}
				}
			}
			if (varToBound != LinearProgram::invalid_variable) {
				/*
				 * Diesen Knoten sofort bearbeiten. So werden schneller obere Schranken gefunden und die LPs können
				 * schneller gelöst werden, da das duale Simplexverfahren genutzt wird und das LP nach dem Hinzufügen
				 * neuer Constraints (bzw. Schranke für Variablen) schnell neu gelöst werden kann
				 * (TODO genauer beschreiben, wenn das in LGO besprochen wurde)
				 */
				branch(varToBound, 1, LinearProgram::lower, node.bounds, fractVal, true);
			}
		}
	}
	currStack.pop_back();
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
 */
void BranchAndCut::branch(variable_id variable, value_t val, LinearProgram::BoundType bound, const SystemBounds& parent,
						  double objValue, bool immediate) {
	assert(parent.getVarCount() == varCount);
	BranchNode node{parent, objValue};
	node.bounds.setBound(variable, bound, val);
	if (immediate || dfs) {
		branchAndBound(node, false);
	} else {
		open.insert(node);
	}
}

/**
 * Erhöht die Zähler für die Constraints, die nicht mit Gleichheit erfüllt sind und setzt die Zähler für die anderen
 * zurück.
 */
void BranchAndCut::countSolutionSlack(const LinearProgram::Solution& sol) {
	const std::vector<double>& slack = sol.getSlack();
	for (size_t constraint = constraintsAtStart; constraint < slack.size(); ++constraint) {
		if (generalTolerance.nonZero(slack[constraint])) {
			++sinceSlack0[constraint - constraintsAtStart];
		} else {
			sinceSlack0[constraint - constraintsAtStart] = 0;
		}
	}
}

/**
 * @param a ein Wert der Zielfunktion
 * @param b ein Wert der Zielfunktion
 * @return true, falls a ein streng besserer Wert der Zielfunktion als b ist
 */
bool BranchAndCut::isBetter(double a, double b) {
	return generalTolerance.less(a, b);
}

/**
 * Setzt die obere Schranke für Branch and Bound
 * @param value Die Werte der Variablen, oder ein leerer Vector
 * @param cost Die Kosten der angegebenen Lösung
 */
void BranchAndCut::setUpperBound(const std::vector<value_t>& value, obj_t cost) {
	upperBound = cost;
	bool variablesValid = static_cast<variable_id>(value.size()) == varCount;
	if (variablesValid) {
		currBest = value;
		std::cout << "Setting upper bound as " << upperBound << std::endl;
	} else {
		std::cout << "Setting upper bound as " << upperBound << " without an example vector" << std::endl;
	}
	if (!open.empty()) {
		size_t removed = 1;
		auto firstRemove = open.end();
		--firstRemove;
		//Letzten Knoten finden, der zu besseren Ergebnissen als der neuen Schranke führen kann
		while (firstRemove != open.begin() && !isBetter(getNextWorse(firstRemove->value, intTolerance), cost)) {
			--firstRemove;
			++removed;
		}
		if (isBetter(getNextWorse(firstRemove->value, intTolerance), cost)) {
			++firstRemove;
			--removed;
		}
		//Alle späteren Knoten entfernen
		open.erase(firstRemove, open.end());
		std::cout << "Removed " << removed << " nodes from the open set" << std::endl;
	}
	if (variablesValid && remover != nullptr) {
		std::vector<variable_id> toRemove = remover->removeOnUpperBound(currBest);
		removeVariables(toRemove);
	}
}

/**
 * Setzt die Variablenbeschränkungen des LPs auf die angegebenen Werte
 */
void BranchAndCut::setupBounds(const SystemBounds& bounds) {
	for (variable_id i = 0; i < varCount; ++i) {
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
	while (i < recentlyRemoved.size()) {
		const LinearProgram::Constraint& constr = recentlyRemoved[i];
		double lhs = constr.evalLHS(sol.getVector());
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

//Entfernt alle Constraints, die seit 10 oder mehr Iterationen nicht mit Gleicheit erfüllt waren, falls es viele gibt
bool BranchAndCut::cleanupOldConstraints() {
	//Es werden erst Constraints entfernt, wenn mehr als doppelt so viele Constraints wie am Anfang der Berechnung existieren
	const double maxRatio = 2;
	//Anzahl an Iterationen, in denen die Constraint nicht mit Gleicheit erfüllt gewesen sein darf, damit sie entfernt wird
	const size_t minSinceEqual = 10;
	//Wenn weniger zu entfernende Constraints als dieser Wert gefunden werden, werden sie noch nicht entfernt
	const size_t minRemove = constraintsAtStart / 10;
	if (problem.getConstraintCount() > constraintsAtStart * maxRatio) {
		std::vector<int> toRemove(problem.getConstraintCount(), 0);
		size_t removeCount = 0;
		std::vector<LinearProgram::Constraint> removed;
		//Alle constraints, die seit 10 oder mehr Iterationen nicht mit Gleichheit erfüllt waren, werden entfernt
		for (size_t i = 0; i < sinceSlack0.size(); ++i) {
			if (sinceSlack0[i] > minSinceEqual) {
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
		std::vector<int> variableMap(varCount, 0);
		for (variable_id r:toRemove) {
			variableMap[r] = 1;
		}
		problem.removeSetVariables(variableMap);
		variable_id newVarCount = varCount - toRemove.size();
		//Variablengrenzen in der offenen Menge updaten
		std::multiset<BranchNode> newOpen;
		for (const BranchNode& old:open) {
			SystemBounds newBounds(old.bounds, variableMap, newVarCount);
			BranchNode newNode{newBounds, old.value};
			//Einen Knoten nicht hinzufügen, wenn er sich nur in den entfernten Variablen von einem anderen Knoten unterschiedet
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
		//Variablengrenzen der aktuellen Knoten updaten
		for (BranchNode *activeNode:currStack) {
			SystemBounds newBounds(activeNode->bounds, variableMap, newVarCount);
			*activeNode = {newBounds, activeNode->value};
		}
		//Einträge zu den entfernten Variablen aus den relevanten vectoren entfernen
		tsp_util::eraseEntries(currentBounds, toRemove);
		tsp_util::eraseEntries(objCoefficients, toRemove);
		tsp_util::eraseEntries(correspondingOrigVar, toRemove);
		for (LinearProgram::Constraint& c:recentlyRemoved) {
			c.deleteVariables(variableMap);
		}
		varCount = newVarCount;
		open = newOpen;
		fractOpt.removeVariables(toRemove);
		assert(varCount == static_cast<variable_id>(currentBounds.size()) &&
			   varCount == static_cast<variable_id>(objCoefficients.size()));
		std::cout << "Removed " << toRemove.size() << " variables, " << varCount << " remaining" << std::endl;
	}
}

value_t& BranchAndCut::VariableBounds::operator[](LinearProgram::BoundType b) {
	if (b == LinearProgram::lower) {
		return min;
	} else {
		return max;
	}
}

bool BranchAndCut::VariableBounds::isFixed() const {
	return max == min;
}

bool BranchAndCut::VariableBounds::operator==(const BranchAndCut::VariableBounds& other) const {
	return min == other.min && max == other.max;
}

bool BranchAndCut::BranchNode::operator<(const BranchAndCut::BranchNode& other) const {
	return value < other.value;
}

bool BranchAndCut::BranchNode::operator==(const BranchAndCut::BranchNode& other) const {
	return bounds == other.bounds;
}

variable_id BranchAndCut::SystemBounds::getVarCount() const {
	return fixLower.size();
}

BranchAndCut::SystemBounds::SystemBounds(variable_id varCount)
		: fixUpper(static_cast<size_t>(varCount), false),
		  fixLower(static_cast<size_t>(varCount), false) {}

/**
 * Erzeugt die SystemBounds, die durch das Entfernen von Variablen aus old entsteht. idMap ist der von LP::removeSet
 * erzeugte vector, der jeder Variablen ihren neuen Index zuordnet
 */
BranchAndCut::SystemBounds::SystemBounds(const BranchAndCut::SystemBounds& old, const std::vector<variable_id>& idMap,
										 variable_id newVarCount) : SystemBounds(newVarCount) {
	assert(idMap.size() == old.fixLower.size());
	for (size_t oldId = 0; oldId < old.fixLower.size(); ++oldId) {
		variable_id newId = idMap[oldId];
		if (newId >= 0) {
			fixUpper[newId] = old.fixUpper[oldId];
			fixLower[newId] = old.fixLower[oldId];
		}
	}
}

BranchAndCut::VariableBounds BranchAndCut::SystemBounds::operator[](variable_id id) const {
	VariableBounds basic{0, 1};
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
	assert(var >= 0 && var < static_cast<variable_id>(fixLower.size()));
	if (b == LinearProgram::lower) {
		fixLower[var] = true;
	} else {
		fixUpper[var] = true;
	}
	assert(!fixUpper[var] || !fixLower[var]);
}

void BranchAndCut::SystemBounds::setBound(variable_id var, LinearProgram::BoundType b, value_t newBound) {
	if (b == LinearProgram::lower) {
		fixUpper[var] = newBound == 1;
	} else {
		fixLower[var] = newBound == 0;
	}
	assert(!fixUpper[var] || !fixLower[var]);
}

bool BranchAndCut::SystemBounds::operator==(const BranchAndCut::SystemBounds& other) const {
	return fixLower == other.fixLower && fixUpper == other.fixUpper;
}