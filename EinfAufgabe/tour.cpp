#include <utility>
//TODO: Warum 2 mal?
#include <utility>
#include "tour.hpp"

Tour::Tour(std::vector<node_id> order, cost_t length, std::string name)
		: length(length), name(std::move(name)), order(std::move(order)) {}

/**
 * Gibt die Tour im TSPLib-Format aus
 */
void Tour::print(std::ostream& out) {
	out << "NAME: " << name << ".tsp.tour\n"
		<< "TYPE: TOUR\n"
		<< "DIMENSION: " << order.size() << "\n"
		<< "TOUR_SECTION\n";
	for (node_id c:order) {
		out << c + 1 << "\n";
	}
	out << "-1\n";
}
