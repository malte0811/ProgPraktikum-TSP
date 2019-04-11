#ifndef COMB_HEURISTIC_HPP
#define COMB_HEURISTIC_HPP

#include <utility>
#include <vector>
#include <tsp_instance.hpp>
#include <contraction_rule.hpp>
#include "blossom_finder.hpp"

class CombHeuristic {
public:
	struct Comb {
		Comb(const Graph& g, const BlossomFinder::Blossom& b,
			 const tsp_util::ContractionMapTSP& contr);

		bool isBlossom() const;

		size_t estimateNonzeroCount() const;

		void invertHandle(city_id cityCount);

		std::vector<city_id> handle;
		std::vector<std::vector<city_id>> teeth;
	};

	explicit CombHeuristic(std::vector<ContractionRule *> rules);

	std::vector<Comb> findViolatedCombs(const TspLpData& lpData, const TSPInstance& inst,
										const std::vector<double>& sol) const;

private:
	std::vector<ContractionRule *> rules;
	lemon::Tolerance<double> tolerance;
};


#endif