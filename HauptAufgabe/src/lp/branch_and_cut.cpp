#include <utility>
#include <branch_and_cut.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <ctime>
#include <queue>
#include <lemon/bin_heap.h>

BranchAndCut::BranchAndCut(LinearProgram& program, const std::vector<CutGenerator*>& gens) :
		problem(program), currBest(program.getVariableCount()), generators(gens),
		fractOpt(program.getVariableCount()), constraintsAtStart(program.getConstraintCount()) {
	if (program.getGoal()==LinearProgram::minimize) {
		upperBound = std::numeric_limits<double>::max();
	} else {
		upperBound = std::numeric_limits<double>::min();
	}
}

void BranchAndCut::solveLP(LinearProgram::Solution& out) {
	bool valid;
	size_t iterations = 0;
	do {
		++iterations;
		if (iterations%64==0) {
			std::cout << "At " << iterations << " iterations, have " << problem.getConstraintCount() << std::endl;
		}
		problem.solve(out);
		if (!out.isValid() || out.getValue()>upperBound) {
			break;
		}
		valid = true;
		for (CutGenerator* gen:generators) {
			if (!gen->validate(problem, out.getVector())) {
				valid = false;
			}
		}
	} while (!valid);
	const double maxRatio = 7;
	if (problem.getConstraintCount()>constraintsAtStart*(maxRatio+1)) {
		std::vector<int> toRemove(problem.getConstraintCount()-constraintsAtStart*maxRatio);
		for (unsigned i = 0; i<toRemove.size(); ++i) {
			toRemove[i] = constraintsAtStart+i;
		}
		problem.removeConstraints(toRemove);
		std::cout << "Removed " << toRemove.size() << " constraints" << std::endl;
	}
}

std::vector<long> BranchAndCut::solve() {
	branchAndBound();
	return currBest;
}

//TODO how much does this actually improve things?
struct BranchInfo {
	int variable;
	double to05;

	bool operator<(const BranchInfo& other) const {
		return to05<other.to05;
	}
};

void BranchAndCut::branchAndBound() {
	solveLP(fractOpt);
	if (!fractOpt.isValid() || !isBetter(fractOpt.getValue(), upperBound)) {
		return;
	}
	std::priority_queue<BranchInfo> possible;
	for (int i = 0; i<problem.getVariableCount(); ++i) {
		long rounded = std::lround(fractOpt[i]);
		double diff = rounded-fractOpt[i];
		if (tolerance.nonZero(diff)) {
			possible.push({i, std::abs(diff-.5)});
		}
	}
	if (possible.empty()) {
		upperBound = fractOpt.getValue();
		for (int i = 0; i<problem.getVariableCount(); ++i) {
			currBest[i] = std::lround(fractOpt[i]);
		}
		//std::cout << "New optimum: " << upperBound << std::endl;
	} else {
		while (!possible.empty() && isBetter(fractOpt.getValue(), upperBound)) {
			int varToBound = possible.top().variable;
			possible.pop();
			long rounded = std::lround(fractOpt[varToBound]);
			double diff = rounded-fractOpt[varToBound];
			long lower;
			if (tolerance.positive(diff)) {
				lower = rounded;
			} else {
				lower = rounded+1;
			}
			long upper = lower-1;
			bound(varToBound, lower, LinearProgram::lower);
			bound(varToBound, upper, LinearProgram::upper);
		}
	}
}

void BranchAndCut::bound(int variable, long val, LinearProgram::BoundType bound) {
	static unsigned calls = 0;
	++calls;
	if (calls%1000==0) {
		std::cout << "Heartbeat: " << std::clock() << std::endl;
		calls = 0;
	}
	double oldBound = problem.getBound(variable, bound);
	problem.setBound(variable, bound, val);
	//std::cout << "Bounding " << variable << " to " << val << " (" << (char) bound << ")" << std::endl;
	branchAndBound();
	//std::cout << "Un-bounding " << variable << "(to " << oldBound << ")" << std::endl;
	problem.setBound(variable, bound, oldBound);
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
