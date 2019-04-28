#ifndef TSP_SOLUTION_HPP
#define TSP_SOLUTION_HPP

#include <tsp_instance.hpp>
#include <lemon/full_graph.h>

class TspLpData;

class TSPSolution {
public:
	TSPSolution() = default;

	TSPSolution(const TSPInstance& inst, const std::vector<bool>& variables, const TspLpData& variableMap);

	TSPSolution(const TSPInstance& inst, const lemon::FullGraph& g, const lemon::FullGraph::EdgeMap<bool>& used);

	TSPSolution(const TSPInstance& inst, std::vector<city_id> order);

	TSPSolution(const TSPInstance& instance, std::istream& input);

	void write(std::ostream& out) const;

	cost_t getCost() const;

	const std::vector<city_id>& getOrder() const;

	TSPSolution opt2() const;

	bool isValid() const;

private:
	void initTourCost();

	void initFromGraph(const lemon::FullGraph& g, const lemon::FullGraph::EdgeMap<bool>& used);

	const TSPInstance *inst = nullptr;
	//Die Reihenfolge der Städte auf der gespeicherten Tour
	std::vector<city_id> order;
	//Die Kosten der Tour
	cost_t cost = 0;
};

#endif