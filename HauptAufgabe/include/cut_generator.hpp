#ifndef CUT_GENERATOR_HPP
#define CUT_GENERATOR_HPP


#include <linear_program.hpp>

class CutGenerator {
public:
	enum CutStatus {
		valid,
		maybe_recalc,
		recalc
	};

	/*
	 * TODO phrasing
	 * Returns true if the given LP solution is valid. Otherwise returns false and adds new constraint that make the
	 * solution invalid.
	 */
	virtual CutStatus validate(LinearProgram& lp, const std::vector<double>& solution) = 0;
};


#endif