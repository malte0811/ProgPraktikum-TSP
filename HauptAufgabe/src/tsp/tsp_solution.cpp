#include <tsp_solution.hpp>
#include <tsp_instance.hpp>
#include <lemon/full_graph.h>
#include <lemon/opt2_tsp.h>
#include <tsp_lp_data.hpp>

/**
 * @param inst Die TSP-Instanz
 * @param variables Die Belegung der LP-Variablen
 */
TSPSolution::TSPSolution(const TSPInstance& inst, const std::vector<bool>& variables, const TspLpData& variableMap)
		: inst(&inst), order(static_cast<size_t>(inst.getCityCount())), variables(variables) {
	city_id previous = 0;
	city_id currentCity = 0;
	size_t indexInTour = 0;
	do {
		if (indexInTour>=order.size()) {
			throw std::runtime_error("Invalid TSP solution, includes a cycle not containing node 1");
		}
		order[indexInTour] = currentCity;
		++indexInTour;
		bool foundNext = false;
		for (city_id next = 0; next<inst.getCityCount(); ++next) {
			variable_id var = variableMap.getVariable(currentCity, next);
			if (next != previous && next != currentCity && var >= 0 && variables[var]) {
				cost += inst.getDistance(currentCity, next);
				previous = currentCity;
				currentCity = next;
				foundNext = true;
				break;
			}
		}
		assert(foundNext);
	} while (currentCity!=0);
	if (indexInTour<order.size()) {
		throw std::runtime_error("Invalid TSP solution, node 1 is in a short cycle");
	}
}

TSPSolution::TSPSolution(const TSPInstance& inst, const std::vector<city_id>& order)
		: inst(&inst), order(order), variables(static_cast<size_t>(inst.getEdgeCount()), false) {
	TspLpData lpData(inst);
	for (city_id i = 0; i < inst.getCityCount(); ++i) {
		variable_id var = lpData.getVariable(order[i], order[(i + 1) % inst.getCityCount()]);
		variables[var] = true;
		cost += lpData.getCost(var);
	}
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
		<< "DIMENSION: " << inst->getCityCount() << "\n"
		<< "TOUR_SECTION\n";
	for (city_id c:getOrder()) {
		out << c+1 << "\n";
	}
	out << "-1\n";
}

TSPSolution TSPSolution::opt2() const {
	using lemon::FullGraph;
	FullGraph g(inst->getCityCount());
	FullGraph::EdgeMap<cost_t> costs(g);
	for (FullGraph::EdgeIt it(g); it != lemon::INVALID; ++it) {
		costs[it] = inst->getDistance(g.id(g.u(it)), g.id(g.v(it)));
	}
	lemon::Opt2Tsp<FullGraph::EdgeMap<cost_t>> opt2(g, costs);
	std::vector<FullGraph::Node> orderGraph(inst->getCityCount());
	for (city_id i = 0; i < inst->getCityCount(); ++i) {
		orderGraph[i] = g(order[i]);
	}
	opt2.run(orderGraph);
	orderGraph = opt2.tourNodes();
	std::vector<city_id> orderInt(inst->getCityCount());
	for (city_id i = 0; i < inst->getCityCount(); ++i) {
		orderInt[i] = g.id(orderGraph[i]);
	}
	return TSPSolution(*inst, orderInt);
}
