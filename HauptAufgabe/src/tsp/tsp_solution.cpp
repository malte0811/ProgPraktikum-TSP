#include <tsp_solution.hpp>
#include <tsp_instance.hpp>

TSPSolution::TSPSolution(const TSPInstance& inst, const std::vector<bool>& variables)
		: inst(inst), variables(variables), order(inst.getSize()) {
	city_id previous = 0;
	city_id currentCity = 0;
	size_t indexInTour = 0;
	//TODO should I handle invalid tours?
	do {
		order[indexInTour] = currentCity;
		++indexInTour;
		for (city_id next = 0; next<inst.getSize(); ++next) {
			if (next!=previous && next!=currentCity && variables[inst.getVariable(currentCity, next)]) {
				cost += inst.getDistance(currentCity, next);
				previous = currentCity;
				currentCity = next;
				break;
			}
		}
	} while (currentCity!=0);
}

cost_t TSPSolution::getCost() const {
	return cost;
}

bool TSPSolution::getVariable(int id) const {
	return variables[id];
}

const std::vector<city_id>& TSPSolution::getOrder() const {
	return order;
}
