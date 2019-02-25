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

	void setUpperBound(const std::vector<long>& value, long cost);

	std::vector<long> solve();

private:

	struct VariableBounds {
		long min, max;

		long& operator[](LinearProgram::BoundType b);
	};

	struct BranchNode {
		std::map<variable_id, VariableBounds> bounds;
		double value;
		LinearProgram::Goal goal;

		bool operator<(const BranchNode& other) const;

		size_t estimateSize() const;
	};

	LinearProgram& problem;
	double upperBound;
	const variable_id varCount;
	const LinearProgram::Goal goal;
	std::vector<long> currBest;
	LinearProgram::Solution fractOpt;
	std::vector<size_t> sinceSlack0;
	std::vector<long> objCoefficients;
	std::set<BranchNode> open;
	std::vector<VariableBounds> defaultBounds;
	std::vector<VariableBounds> currentBounds;
	size_t openSize = 0;
	const std::vector<CutGenerator*> generators;
	const size_t constraintsAtStart;
	const lemon::Tolerance<double> generalTolerance;
	const lemon::Tolerance<double> intTolerance;

	static bool isBetter(double a, double b, LinearProgram::Goal goal);

	void branchAndBound(BranchNode& node, bool dfs);

	void branch(int variable, long val, LinearProgram::BoundType bound,
				const std::map<variable_id, VariableBounds>& parent, double objValue, bool immediate, bool dfs);

	void solveLP(LinearProgram::Solution& out);

	void countSolutionSlack(const LinearProgram::Solution& sol);

	void setupBounds(std::map<variable_id, VariableBounds> bounds);
};


#endif