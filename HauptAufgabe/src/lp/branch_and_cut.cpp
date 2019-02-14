#include <utility>
#include <branch_and_cut.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <ctime>
#include <queue>

BranchAndCut::BranchAndCut(LinearProgram& program, const std::vector<CutGenerator*>& gens) :
		problem(program), currBest(static_cast<size_t>(program.getVariableCount())),
		fractOpt(static_cast<size_t>(program.getVariableCount())),
		generators(gens),
		constraintsAtStart(static_cast<size_t>(program.getConstraintCount())) {
	if (program.getGoal()==LinearProgram::minimize) {
		upperBound = std::numeric_limits<double>::max();
	} else {
		upperBound = -std::numeric_limits<double>::max();
	}
}

void BranchAndCut::solveLP(LinearProgram::Solution& out) {
	CutGenerator::CutStatus solutionValid;
	double oldVal = 0;
	size_t slowIterations = 0;
	do {
		problem.solve(out);
		countSolutionSlack();
		if (!out.isValid() || !isBetter(out.getValue(), upperBound)) {
			break;
		}
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
		sinceSlack0.resize(problem.getConstraintCount()-constraintsAtStart, 0);
		if (slowIterations>6 && solutionValid==CutGenerator::maybe_recalc) {
			solutionValid = CutGenerator::valid;
		}
	} while (solutionValid!=CutGenerator::valid);
	const double maxRatio = 3;
	if (problem.getConstraintCount()>constraintsAtStart*maxRatio) {
		std::vector<int> toRemove;
		for (size_t i = 0; i<sinceSlack0.size(); ++i) {
			if (sinceSlack0[i]>10) {
				toRemove.push_back(i+constraintsAtStart);
			}
		}
		if (!toRemove.empty()) {
			std::vector<size_t> since0New;
			since0New.reserve(sinceSlack0.size()-toRemove.size());
			size_t nextInRemove = 0;
			for (size_t i = 0; i<sinceSlack0.size(); ++i) {
				if (nextInRemove<toRemove.size() && toRemove[nextInRemove]==(i+constraintsAtStart)) {
					++nextInRemove;
				} else {
					since0New.push_back(sinceSlack0[i]);
				}
			}
			sinceSlack0 = std::move(since0New);
			problem.removeConstraints(toRemove);
			//std::cout << "Removed " << toRemove.size() << " constraints" << std::endl;
		}
	}
}

static size_t boundSteps = 0;

std::vector<long> BranchAndCut::solve() {
	boundSteps = 0;
	branchAndBound();
	std::cout << boundSteps << std::endl;
	return currBest;
}

void BranchAndCut::branchAndBound() {
	solveLP(fractOpt);
	if (!fractOpt.isValid() || !isBetter(fractOpt.getValue(), upperBound)) {
		return;
	}
	variable_id varToBound = -1;
	double optDist = 1;
	for (variable_id i = 0; i<problem.getVariableCount(); ++i) {
		long rounded = std::lround(fractOpt[i]);
		double diff = rounded-fractOpt[i];
		if (tolerance.nonZero(diff)) {
			double dist05 = std::abs(std::abs(diff)-.5);
			if (dist05<optDist) {
				optDist = dist05;
				varToBound = i;
			}
		}
	}
	if (varToBound<0) {
		upperBound = fractOpt.getValue();
		for (int i = 0; i<problem.getVariableCount(); ++i) {
			currBest[i] = std::lround(fractOpt[i]);
		}
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
		bound(varToBound, lower, LinearProgram::lower);
		bound(varToBound, upper, LinearProgram::upper);
	}
}

void BranchAndCut::bound(int variable, long val, LinearProgram::BoundType bound) {
	++boundSteps;
	double oldBound = problem.getBound(variable, bound);
	problem.setBound(variable, bound, val);
	branchAndBound();
	problem.setBound(variable, bound, oldBound);
}

void BranchAndCut::countSolutionSlack() {
	std::vector<double> slack = problem.getSlack();
	for (size_t constraint = constraintsAtStart; constraint<slack.size(); ++constraint) {
		if (tolerance.nonZero(slack[constraint])) {
			++sinceSlack0[constraint-constraintsAtStart];
		} else {
			sinceSlack0[constraint-constraintsAtStart] = 0;
		}
	}
}

bool BranchAndCut::isBetter(double a, double b) {
	if (problem.getGoal()==LinearProgram::maximize) {
		return std::round(b)<std::round(a);
	} else {
		return std::round(a)<std::round(b);
	}
}

void BranchAndCut::setUpperBound(const std::vector<long>& value, double cost) {
	upperBound = cost;
	currBest = value;
}
