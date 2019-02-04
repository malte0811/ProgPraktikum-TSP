#ifndef HAUPTAUFGABE_BRANCH_AND_CUT_HPP
#define HAUPTAUFGABE_BRANCH_AND_CUT_HPP


#include <lemon/tolerance.h>
#include <linear_program.hpp>
#include <cut_generator.hpp>

class BranchAndCut {
public:
	BranchAndCut(LinearProgram& program, const std::vector<CutGenerator*>& gens);

	void setUpperBound(const std::vector<long>& value, double cost);

	std::vector<long> solve();

private:
	LinearProgram& problem;
	double upperBound;
	std::vector<long> currBest;
	LinearProgram::Solution fractOpt;
	const std::vector<CutGenerator*> generators;
	lemon::Tolerance<double> tolerance;

	bool isBetter(double a, double b);

	void branchAndBound();

	void bound(int variable, long val, LinearProgram::BoundType bound);

	void solveLP(LinearProgram::Solution& out);
};


#endif //HAUPTAUFGABE_BRANCH_AND_CUT_HPP
