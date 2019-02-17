#ifndef BRANCH_AND_CUT_HPP
#define BRANCH_AND_CUT_HPP


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
	std::vector<size_t> sinceSlack0;
	std::vector<long> objCoefficients;
	const std::vector<CutGenerator*> generators;
	const size_t constraintsAtStart;
	const lemon::Tolerance<double> tolerance;

	bool isBetter(double a, double b);

	void branchAndBound();

	void bound(int variable, long val, LinearProgram::BoundType bound);

	void solveLP(LinearProgram::Solution& out);

	void countSolutionSlack();
};


#endif