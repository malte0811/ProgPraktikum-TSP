#ifndef BRANCH_AND_CUT_HPP
#define BRANCH_AND_CUT_HPP


#include <lemon/tolerance.h>
#include <linear_program.hpp>
#include <cut_generator.hpp>
#include <queue>
#include <stack>
#include <iostream>
#include <array>
#include <map>
#include <set>

class BranchAndCut {
public:
	BranchAndCut(LinearProgram& program, const std::vector<CutGenerator*>& gens);

	void setUpperBound(const std::vector<long>& value, double cost);

	std::vector<long> solve();

private:

	using VariableBounds = std::array<long, 2>;

	struct BranchNode {
		std::map<variable_id, VariableBounds> bounds;
		double value;
		size_t level;
		LinearProgram::Goal goal;

		bool operator<(const BranchNode& other) const {
			return value<other.value;
		}
	};

	LinearProgram& problem;
	double upperBound;
	std::vector<long> currBest;
	LinearProgram::Solution fractOpt;
	std::vector<size_t> sinceSlack0;
	std::vector<long> objCoefficients;
	std::set<BranchNode> open;
	std::vector<VariableBounds> defaultBounds;
	const std::vector<CutGenerator*> generators;
	const size_t constraintsAtStart;
	const lemon::Tolerance<double> tolerance;

	static bool isBetter(double a, double b, LinearProgram::Goal goal);

	void branchAndBound(const BranchNode& node, bool setup);

	void bound(int variable, long val, LinearProgram::BoundType bound,
			   const std::map<variable_id, VariableBounds>& parent,
			   size_t level, double objValue, bool immediate);

	void solveLP(LinearProgram::Solution& out);

	void countSolutionSlack();

	void setupBounds(std::map<variable_id, VariableBounds> bounds);
};


#endif