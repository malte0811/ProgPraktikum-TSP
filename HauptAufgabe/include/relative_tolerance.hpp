#ifndef RELATIVE_TOLERANCE_HPP
#define RELATIVE_TOLERANCE_HPP


#include <cassert>

class RelativeTolerance {
public:
	RelativeTolerance() = default;

	static const double eps;

	bool less(double a, double b) const;

	bool different(double a, double b) const;

	bool positive(double a) const;

	bool negative(double a) const;

	bool nonZero(double a) const;
};

#endif
