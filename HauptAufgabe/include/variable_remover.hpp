#ifndef VARIABLE_REMOVER_HPP
#define VARIABLE_REMOVER_HPP

#include <vector>
#include <linear_program.hpp>
#include <branch_and_cut.hpp>

class VariableRemover {
public:
	virtual std::vector<variable_id> removeOnUpperBound(const std::vector<value_t>& variables) = 0;

	virtual std::vector<variable_id> removeOnRootSolution(const LinearProgram::Solution& rootSol) = 0;
};

#endif