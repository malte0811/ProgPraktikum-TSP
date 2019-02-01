#include <utility>
#include <branch_and_cut.hpp>
#include <cmath>
#include <iostream>
#include <limits>

BranchAndCut::BranchAndCut(LinearProgram& p, double initialUpperBound, std::vector<long> initialOpt,
						   const std::vector<CutGenerator*>& gens)
		: problem(p), upperBound(initialUpperBound), currBest(std::move(initialOpt)), generators(gens) {}

BranchAndCut::BranchAndCut(LinearProgram& program, const std::vector<CutGenerator*>& gens) :
		BranchAndCut(program, 0, std::vector<long>(static_cast<size_t>(program.getVariableCount())), gens) {
	if (program.getGoal()==LinearProgram::minimize) {
		upperBound = std::numeric_limits<double>::max();
	} else {
		upperBound = std::numeric_limits<double>::min();
	}
}

std::vector<double> BranchAndCut::solveLP() {
	std::vector<double> ret;
	bool valid;
	do {
		ret = problem.solve();
		if (ret.empty() || ret.back()>upperBound) {
			break;
		}
		valid = true;
		for (CutGenerator* gen:generators) {
			if (!gen->validate(problem, ret)) {
				valid = false;
				break;
			}
		}
	} while (!valid);
	return ret;
}

std::vector<long> BranchAndCut::solve() {
	branchAndBound();
	return currBest;
}

void BranchAndCut::branchAndBound() {
	std::vector<double> opt = solveLP();
	if (opt.empty() || !isBetter(opt.back(), upperBound)) {
		return;
	}
	bool integer = true;
	for (int i = 0; i<problem.getVariableCount(); ++i) {
		//TODO heuristics for choosing which variable to bound?
		long rounded = std::lround(opt[i]);
		double diff = rounded-opt[i];
		if (diff>eps) {
			bound(i, rounded, LinearProgram::lower);
			bound(i, rounded-1, LinearProgram::upper);
			integer = false;
		} else if (diff<-eps) {
			bound(i, rounded+1, LinearProgram::lower);
			bound(i, rounded, LinearProgram::upper);
			integer = false;
		}
	}
	if (integer) {
		upperBound = opt.back();
		for (int i = 0; i<problem.getVariableCount(); ++i) {
			currBest[i] = std::lround(opt[i]);
		}
	}
}

void BranchAndCut::bound(int variable, long val, LinearProgram::BoundType bound) {
	double oldBound = problem.getBound(variable, bound);
	problem.setBound(variable, bound, val);
	std::cout << "Bounding " << variable << " to " << val << " (" << (char) bound << ")" << std::endl;
	branchAndBound();
	std::cout << "Un-bounding " << variable << "(to " << oldBound << ")" << std::endl;
	problem.setBound(variable, bound, oldBound);
}

bool BranchAndCut::isBetter(double a, double b) {
	if (problem.getGoal()==LinearProgram::maximize) {
		return a>b;
	} else {
		return a<b;
	}
}
