#ifndef VARIABLE_REMOVER_HPP
#define VARIABLE_REMOVER_HPP

#include <vector>
#include <linear_program.hpp>
#include <branch_and_cut.hpp>

class VariableRemover {
public:
	virtual std::vector<variable_id> removeVariables(long newBound, const std::vector<long>& variables) = 0;
};

#endif