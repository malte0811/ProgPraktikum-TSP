#include <utility>
#include <branch_and_cut.hpp>
#include <cmath>
#include <iostream>
#include <limits>

BranchAndCut::BranchAndCut(LinearProgram& p, double initialUpperBound, std::vector<long> initialOpt,
						   const std::vector<CutGenerator*>& gens)
		: problem(p), upperBound(initialUpperBound), currBest(std::move(initialOpt)), generators(gens),
		  fractOpt(p.getVariableCount()) {}

BranchAndCut::BranchAndCut(LinearProgram& program, const std::vector<CutGenerator*>& gens) :
		BranchAndCut(program, 0, std::vector<long>(static_cast<size_t>(program.getVariableCount())), gens) {
	if (program.getGoal()==LinearProgram::minimize) {
		upperBound = std::numeric_limits<double>::max();
	} else {
		upperBound = std::numeric_limits<double>::min();
	}
}

void BranchAndCut::solveLP(LinearProgram::Solution& out) {
	bool valid;
	do {
		problem.solve(out);
		if (!out.isValid() || out.getValue()>upperBound) {
			break;
		}
		valid = true;
		/*
		int hash = 0;
		for (double d:out.getVector()) {
			hash *= 37;
			hash += std::lround(d);
		}
		std::cout << hash << std::endl;
		*/
		for (CutGenerator* gen:generators) {
			if (!gen->validate(problem, out.getVector())) {
				valid = false;
				break;
			}
		}
	} while (!valid);
}

std::vector<long> BranchAndCut::solve() {
	branchAndBound();
	return currBest;
}

void BranchAndCut::branchAndBound() {
	solveLP(fractOpt);
	if (!fractOpt.isValid() || !isBetter(fractOpt.getValue(), upperBound)) {
		return;
	}
	bool integer = true;
	for (int i = 0; i<problem.getVariableCount(); ++i) {
		//TODO heuristics for choosing which variable to bound?
		long rounded = std::lround(fractOpt[i]);
		double diff = rounded-fractOpt[i];
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
		upperBound = fractOpt.getValue();
		for (int i = 0; i<problem.getVariableCount(); ++i) {
			currBest[i] = std::lround(fractOpt[i]);
		}
		std::cout << "New optimum: " << upperBound << std::endl;
	}
}

void BranchAndCut::bound(int variable, long val, LinearProgram::BoundType bound) {
	double oldBound = problem.getBound(variable, bound);
	problem.setBound(variable, bound, val);
	//std::cout << "Bounding " << variable << " to " << val << " (" << (char) bound << ")" << std::endl;
	branchAndBound();
	//std::cout << "Un-bounding " << variable << "(to " << oldBound << ")" << std::endl;
	problem.setBound(variable, bound, oldBound);
}

bool BranchAndCut::isBetter(double a, double b) {
	if (problem.getGoal()==LinearProgram::maximize) {
		return a>b;
	} else {
		return a<b;
	}
}
