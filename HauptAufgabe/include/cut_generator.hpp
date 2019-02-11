#ifndef CUT_GENERATOR_HPP
#define CUT_GENERATOR_HPP


#include <linear_program.hpp>

class CutGenerator {
public:
	/*
	 * TODO phrasing
	 * Returns true if the given LP solution is valid. Otherwise returns false and adds new constraint that make the
	 * solution invalid.
	 */
	virtual bool validate(LinearProgram& lp, const std::vector<double>& solution) = 0;
};


#endif