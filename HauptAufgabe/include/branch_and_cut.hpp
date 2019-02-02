#ifndef HAUPTAUFGABE_BRANCH_AND_CUT_HPP
#define HAUPTAUFGABE_BRANCH_AND_CUT_HPP


#include "linear_program.hpp"
#include "cut_generator.hpp"

class BranchAndCut {
public:
	BranchAndCut(LinearProgram& p, double initialUpperBound, std::vector<long> initialOpt,
				 const std::vector<CutGenerator*>& gens);

	BranchAndCut(LinearProgram& program, const std::vector<CutGenerator*>& gens);

	std::vector<long> solve();

private:
	LinearProgram& problem;
	double upperBound;
	std::vector<long> currBest;
	LinearProgram::Solution fractOpt;
	const std::vector<CutGenerator*> generators;
	static constexpr double eps = 0.01;

	bool isBetter(double a, double b);

	void branchAndBound();

	void bound(int variable, long val, LinearProgram::BoundType bound);

	void solveLP(LinearProgram::Solution& out);
};


#endif //HAUPTAUFGABE_BRANCH_AND_CUT_HPP
