#ifndef BRANCH_AND_CUT_HPP
#define BRANCH_AND_CUT_HPP


#include <lemon/tolerance.h>
#include <linear_program.hpp>
#include <cut_generator.hpp>
#include <queue>
#include <stack>
#include <iostream>

class BranchAndCut {
public:
	BranchAndCut(LinearProgram& program, const std::vector<CutGenerator*>& gens);

	void setUpperBound(const std::vector<long>& value, double cost);

	std::vector<long> solve();

private:

	struct Bound {
		variable_id var;
		long value;
		LinearProgram::BoundType type;
	};

	struct BranchNode {
		std::vector<Bound> bounds;
		double value;
		size_t level;
		LinearProgram::Goal goal;
	};

	LinearProgram& problem;
	double upperBound;
	std::vector<long> currBest;
	LinearProgram::Solution fractOpt;
	std::vector<size_t> sinceSlack0;
	std::vector<long> objCoefficients;
	std::stack<BranchNode> open;
	const std::vector<CutGenerator*> generators;
	const size_t constraintsAtStart;
	const lemon::Tolerance<double> tolerance;

	static bool isBetter(double a, double b, LinearProgram::Goal goal);

	void branchAndBound(const BranchNode& node, bool setup);

	void bound(int variable, long val, LinearProgram::BoundType bound, const std::vector<Bound>& parent,
			   size_t level, double objValue, bool immediate);

	void solveLP(LinearProgram::Solution& out);

	void countSolutionSlack();

	std::vector<double> setupBounds(const std::vector<Bound>& bounds);

	void cleanupBounds(const std::vector<Bound>& bounds, const std::vector<double>& oldBounds);
};


#endif