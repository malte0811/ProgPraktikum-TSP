#include <tsp_solution.hpp>
#include <tsp_instance.hpp>

TSPSolution::TSPSolution(const TSPInstance& inst, const std::vector<bool>& variables)
		: inst(&inst), order(static_cast<size_t>(inst.getSize())), variables(variables) {
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

/**
 * @return Die Kosten der Lösung
 */
cost_t TSPSolution::getCost() const {
	return cost;
}

/**
 * @param id Die ID einer Variablen des LP's
 * @return true, falls die Variable den Wert 1 hat
 */
bool TSPSolution::getVariable(int id) const {
	return variables[id];
}

/**
 * @return Die Städte in der Reihenfolge, in der sie in dieser Lösung besucht werden
 */
const std::vector<city_id>& TSPSolution::getOrder() const {
	return order;
}

/**
 * Gibt die Lösung im TSPLIB-Format in out aus
 */
void TSPSolution::write(std::ostream& out) const {
	if (inst==nullptr) {
		throw std::runtime_error("Tried to write invalid TSP solution!");
	}
	out << "NAME: " << inst->getName() << ".tsp.tour\n"
		<< "TYPE: TOUR\n"
		<< "DIMENSION: " << inst->getSize() << "\n"
		<< "TOUR_SECTION\n";
	for (city_id c:getOrder()) {
		out << c+1 << "\n";
	}
	out << "-1\n";
}