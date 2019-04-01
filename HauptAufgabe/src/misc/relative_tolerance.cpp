#include <relative_tolerance.hpp>
#include <lemon/tolerance.h>

RelativeTolerance::RelativeTolerance() : RelativeTolerance(lemon::Tolerance<double>::defaultEpsilon()) {}

RelativeTolerance::RelativeTolerance(double eps) : eps(eps) {}

bool RelativeTolerance::less(double a, double b) const {
	if (a>=b) {
		return false;
	}
	if (a<0) {
		if (b<0) {
			return a<b*(1+eps);
		} else {
			return a+eps<b;
		}
	} else {
		return a*(1+eps)<b;
	}
}

bool RelativeTolerance::different(double a, double b) const {
	return less(a, b) || less(b, a);
}

bool RelativeTolerance::positive(double a) const {
	return eps<a;
}

bool RelativeTolerance::negative(double a) const {
	return a<-eps;
}

bool RelativeTolerance::nonZero(double a) const {
	return negative(a) || positive(a);
}
