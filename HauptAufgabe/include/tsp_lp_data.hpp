#ifndef TSP_LP_DATA_HPP
#define TSP_LP_DATA_HPP


#include <lemon/bfs.h>
#include <variable_remover.hpp>
#include <tsp_instance.hpp>
#include <one_tree.hpp>
#include "tsp_solution.hpp"

class TspLpData : public VariableRemover {
public:
	using Edge = std::pair<city_id, city_id>;

	explicit TspLpData(const TSPInstance& inst);

	std::vector<variable_id> removeVariables(const std::vector<value_t>& variables) override;

	std::vector<variable_id> removeVariables(const TSPSolution& solution);

	inline cost_t getCost(variable_id var) const;

	inline Edge getEdge(variable_id var) const;

	void setupBasicLP(LinearProgram& lp) const;

	inline variable_id getVariable(city_id a, city_id b) const;

	const TSPSolution& getUpperBound() const;

	void setupLowerBounds();

	inline variable_id getVariableCount() const;
private:
	const TSPInstance& inst;
	std::vector<std::pair<city_id, city_id>> variableToEdge;
	std::vector<std::vector<variable_id>> edgeToVariable;
	std::vector<double> removalBound;
	std::vector<bool> shouldHaveRemoved;
	TSPSolution upperBound;
};

TspLpData::Edge TspLpData::getEdge(variable_id var) const {
	return variableToEdge[var];
}

cost_t TspLpData::getCost(variable_id var) const {
	Edge e = getEdge(var);
	return inst.getDistance(e.first, e.second);
}

variable_id TspLpData::getVariable(city_id a, city_id b) const {
	//TODO assert(a!=b);
	if (a == b) return -1;
	if (a < b) {
		std::swap(a, b);
	}
	return edgeToVariable[a - 1][b];
}

variable_id TspLpData::getVariableCount() const {
	return variableToEdge.size();
}

#endif