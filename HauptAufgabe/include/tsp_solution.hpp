#ifndef TSP_SOLUTION_HPP
#define TSP_SOLUTION_HPP

#include <tsp_instance.hpp>

class TspLpData;

class TSPSolution {
public:
	TSPSolution() = default;

	TSPSolution(const TSPInstance& inst, const std::vector<bool>& variables, const TspLpData& variableMap);

	TSPSolution(const TSPInstance& inst, const std::vector<city_id>& order);

	void write(std::ostream& out) const;

	cost_t getCost() const;

	bool getVariable(variable_id id) const;

	const std::vector<city_id>& getOrder() const;

	TSPSolution opt2() const;

private:
	const TSPInstance* inst = nullptr;
	//Die Reihenfolge der Städte auf der gespeicherten Tour
	std::vector<city_id> order;
	//Die Belegung der LP-Variablen
	std::vector<bool> variables;
	//Die Kosten der Tour
	cost_t cost = 0;
};

#endif