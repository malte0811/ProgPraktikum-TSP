#ifndef COMB_HEURISTIC_HPP
#define COMB_HEURISTIC_HPP

#include <utility>
#include <vector>
#include <tsp_instance.hpp>
#include <contraction_rule.hpp>

class CombHeuristic {
public:
	struct Comb {
		std::vector<city_id> handle;
		std::vector<std::vector<city_id>> teeth;
	};

	explicit CombHeuristic(std::vector<ContractionRule *> rules);

	std::vector<Comb> findViolatedCombs(const TspLpData& lpData, const TSPInstance& inst,
										const std::vector<double>& sol);

private:
	std::vector<ContractionRule *> rules;
	lemon::Tolerance<double> tolerance;
};


#endif