#ifndef TWO_MATCHING_CUT_GEN_HPP
#define TWO_MATCHING_CUT_GEN_HPP

#include <lemon/tolerance.h>
#include <cut_generator.hpp>
#include <vector>

class LinearProgram;

class TSPInstance;

class TspLpData;

class TwoMatchingCutGen : public CutGenerator {
public:
	explicit TwoMatchingCutGen(const TSPInstance& inst, const TspLpData& lpData, bool contract);

	CutStatus validate(LinearProgram& lp, const std::vector<double>& solution, CutStatus currentStatus) override;

private:
	const TSPInstance& tsp;

	const TspLpData& lpData;

	const bool enableContraction;

	const lemon::Tolerance<double> tolerance;
};


#endif
