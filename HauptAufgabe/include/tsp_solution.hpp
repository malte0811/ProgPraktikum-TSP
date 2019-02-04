#ifndef TSP_SOLUTION_HPP
#define TSP_SOLUTION_HPP


#include "tsp_instance.hpp"

class TSPSolution {
public:
	TSPSolution(const TSPInstance& inst, const std::vector<bool>& variables);

	//TODO do I need this one anywhere?
	TSPSolution(const TSPInstance& inst, std::vector<city_id> order);

	void write(std::ostream& out) const;

	cost_t getCost() const;

	bool getVariable(variable_id id) const;

	const std::vector<city_id>& getOrder() const;

private:
	const TSPInstance& inst;
	std::vector<city_id> order;
	std::vector<bool> variables;
	cost_t cost = 0;
};

#endif