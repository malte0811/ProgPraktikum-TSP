#ifndef COMB_HEURISTIC_HPP
#define COMB_HEURISTIC_HPP

#include <utility>
#include <vector>
#include <tsp_instance.hpp>
#include <contraction_rule.hpp>
#include "blossom_finder.hpp"

//Eine Instanz der Generic Comb Heruistic aus Grötschel, Holland 1991
class CombHeuristic {
public:
	struct Comb {
		Comb(const Graph& g, const BlossomFinder::Blossom& b,
			 const tsp_util::ContractionMapTSP& contr);

		//Wahr, falls der Kamm eine Blüte ist, d.h. alle Zinken Größe 2 haben
		bool isBlossom() const;

		//Vereinfacht den Kamm, z.B. durch das Entfernen langer Ketten von 1-Kanten aus dem Griff
		void simplify(const TspLpData& lpData, const std::vector<double>& solution);

		//Prüft, dass der Kamm "gültig" (d.h. tatsächlich ein Kamm) ist
		void validate(city_id cityCount) const;

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