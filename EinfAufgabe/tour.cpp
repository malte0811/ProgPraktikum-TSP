#include <utility>
#include "tour.hpp"

Tour::Tour(std::vector<node_id> order, cost_t length, const std::string& name)
		: order(std::move(order)), name(name), length(length) {}

/**
 * Gibt die Tour im TSPLib-Format aus
 */
void Tour::print(std::ostream& out) {
	out << "NAME: " << name << ".tsp.tour\n"
		<< "TYPE: TOUR\n"
		<< "DIMENSION: " << order.size() << "\n"
		<< "TOUR_SECTION\n";
	for (node_id c:order) {
		out << c+1 << "\n";
	}
	out << "-1\n";
}
