#include <utility>
#include <branch_and_cut.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <ctime>
#include <queue>

BranchAndCut::BranchAndCut(LinearProgram& program, const std::vector<CutGenerator*>& gens, size_t maxOpenSize) :
		problem(program), varCount(program.getVariableCount()), goal(program.getGoal()),
		currBest(static_cast<size_t>(varCount)), fractOpt(static_cast<size_t>(varCount)),
		objCoefficients(static_cast<size_t>(varCount)), maxOpenSize(maxOpenSize), generators(gens),
		constraintsAtStart(static_cast<size_t>(program.getConstraintCount())),
		defaultBounds(static_cast<size_t>(varCount)),
		currentBounds(static_cast<size_t>(varCount)),
		intTolerance(0.01) {
	if (goal==LinearProgram::minimize) {
		upperBound = std::numeric_limits<double>::max();
	} else {
		upperBound = -std::numeric_limits<double>::max();
	}
	std::vector<double> obj = program.getObjective();
	for (size_t i = 0; i<obj.size(); ++i) {
		objCoefficients[i] = std::lround(obj[i]);
	}
	for (variable_id i = 0; i<varCount; ++i) {
		for (LinearProgram::BoundType type:{LinearProgram::lower, LinearProgram::upper}) {
			defaultBounds[i][type] = std::lround(problem.getBound(i, type));
		}
		currentBounds[i] = defaultBounds[i];
	}
}

double tolerantCeil(double in, lemon::Tolerance<double> tol) {
	double ceil = std::ceil(in);
	if (tol.different(ceil-in, 1)) {
		return ceil;
	} else {
		return ceil-1;
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
 * Amd Ende werden Constraints entfernt, die seit mindestens 10 Iterationen nicht mehr mit Gleichheit erfüllt waren.
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
		 * TODO Annahme: Koeffs der Zielfunktion sind ganzzahlig
		 * TODO ceil durch floor ersetzen, falls max. wird
		 */
		if (!out.isValid() || !isBetter(tolerantCeil(out.getValue(), intTolerance), upperBound, goal)) {
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
			//std::cout << "Removing " << toRemove.size() << std::endl;
		}
	}
}

/**
 * @return eine optimale ganzzahlige Lösung des LP, die von alle Cut-Generatoren akzeptiert wird.
 */
std::vector<long> BranchAndCut::solve() {
	BranchNode initNode{SystemBounds(this), 0, goal};
	open.insert(initNode);
	openSize += initNode.estimateSize();
	while (!open.empty()) {
		auto it = open.begin();
		BranchNode next = *it;
		open.erase(it);
		openSize -= next.estimateSize();
		branchAndBound(next, false);
		std::cout << openSize << ", size: " << open.size() << std::endl;
		if (maxOpenSize>0) {
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
	if (fractOpt.isValid() && isBetter(tolerantCeil(fractOpt.getValue(), intTolerance), upperBound, goal)) {
		variable_id varToBound = -1;
		double optDist = 1;
		long varWeight = 0;
		const std::vector<double>& reducedCosts = fractOpt.getReducedCosts();
		size_t nonIntCount = 0;
		for (variable_id i = 0; i<varCount; ++i) {
			long rounded = std::lround(fractOpt[i]);
			double diff = rounded-fractOpt[i];
			if (intTolerance.nonZero(diff)) {
				double dist05 = std::abs(std::abs(diff)-.5);
				long cost = objCoefficients[i];
				if (generalTolerance.less(dist05, optDist) ||
					(!generalTolerance.less(optDist, dist05) && cost>varWeight)) {
					optDist = dist05;
					varToBound = i;
					varWeight = cost;
				}
				++nonIntCount;
			} else {
				LinearProgram::BoundType upper = LinearProgram::upper, lower = LinearProgram::lower;
				if (goal==LinearProgram::maximize) {
					std::swap(upper, lower);
				}
				if (node.bounds[i][upper]!=node.bounds[i][lower]) {
					if (rounded==node.bounds[i][upper] && reducedCosts[i]<-(upperBound-fractOpt.getValue())) {
						problem.setBound(i, lower, rounded);
						node.bounds.fix(i, upper);
						currentBounds[i] = node.bounds[i];
					} else if (rounded==node.bounds[i][lower] && reducedCosts[i]>upperBound-fractOpt.getValue()) {
						problem.setBound(i, upper, rounded);
						node.bounds.fix(i, lower);
						currentBounds[i] = node.bounds[i];
					}
				}
			}
		}
		if (varToBound<0) {
			std::vector<long> newOpt(varCount);
			for (int i = 0; i<varCount; ++i) {
				newOpt[i] = std::lround(fractOpt[i]);
			}
			setUpperBound(newOpt, std::lround(fractOpt.getValue()));
		} else {
			long rounded = std::lround(fractOpt[varToBound]);
			double diff = rounded-fractOpt[varToBound];
			long lower;
			if (diff>0) {
				lower = rounded;
			} else {
				lower = rounded+1;
			}
			long upper = lower-1;
			double fractVal = fractOpt.getValue();
			branch(varToBound, lower, LinearProgram::lower, node.bounds, fractVal,
				   dfs || openSize<maxOpenSize, dfs);
			branch(varToBound, upper, LinearProgram::upper, node.bounds, fractVal,
				   dfs || (nonIntCount<10 && openSize<maxOpenSize), dfs);
		}
	}
}

/**
 * Beschränkt den zulässgen Bereich einer Variablen wie durch die Parameter beschrieben, ruft dann branchAndBound auf
 * und setzt den zulässigen Bereich wieder auf den ursprünglichen Wert.
 * @param variable Die zu beschränkende Variable
 * @param val Der neue Wert der Beschränkung
 * @param bound Die "Richtung", in der die variable beschränkt werden soll
 * TODO update comment
 */
void BranchAndCut::branch(int variable, long val, LinearProgram::BoundType bound,
						  const SystemBounds& parent, double objValue, bool immediate, bool dfs) {
	BranchNode node{parent, objValue, goal};
	if (!node.bounds.bounds.count(variable)) {
		node.bounds.bounds[variable] = defaultBounds[variable];
	}
	node.bounds.bounds[variable][bound] = val;
	if (immediate) {
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
	auto it = open.cbegin();
	while (it!=open.cend() && !isBetter(tolerantCeil(it->value, intTolerance), cost, goal)) {
		openSize -= it->estimateSize();
		it = open.erase(it);
	}
}

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

size_t BranchAndCut::BranchNode::estimateSize() const {
	return sizeof(*this)+bounds.bounds.size()*sizeof(*bounds.bounds.begin())+
		   2*bounds.fixLower.size()/CHAR_BIT;
}

bool BranchAndCut::BranchNode::operator<(const BranchAndCut::BranchNode& other) const {
	return value<other.value;
}

BranchAndCut::SystemBounds::SystemBounds(BranchAndCut* owner) : owner(owner),
																fixLower(owner->varCount, false),
																fixUpper(owner->varCount, false) {}

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

void BranchAndCut::SystemBounds::fix(variable_id var, LinearProgram::BoundType b) {
	if (b==LinearProgram::lower) {
		fixLower[var] = true;
	} else {
		fixUpper[var] = true;
	}
}
