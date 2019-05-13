#include <lemon/core.h>
#include <lemon/full_graph.h>
#include <lemon/opt2_tsp.h>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tsp_instance.hpp>
#include <tsp_lp_data.hpp>
#include <tsp_solution.hpp>
#include <tsp_utils.hpp>
#include <utility>
#include <vector>
#include <linear_program.hpp>

/**
 * @param inst Die TSP-Instanz
 * @param variables Die Belegung der LP-Variablen
 * @param variableMap gibt an, welche LP-Variable welcher Kante entspricht
 */
TSPSolution::TSPSolution(const TSPInstance& inst, const std::vector<bool>& variables, const TspLpData& variableMap)
		: inst(&inst), order(static_cast<size_t>(inst.getCityCount())) {
	lemon::FullGraph g(inst.getCityCount());
	lemon::FullGraph::EdgeMap<bool> used(g);
	for (variable_id i = 0; i < variableMap.getVariableCount(); ++i) {
		if (variables[i]) {
			const TspLpData::Edge& e = variableMap.getEdge(i);
			lemon::FullGraph::Node endA = g(e.first);
			lemon::FullGraph::Node endB = g(e.second);
			used[g.edge(endA, endB)] = true;
		}
	}
	initFromGraph(g, used);
}

void TSPSolution::initFromGraph(const lemon::FullGraph& g, const lemon::FullGraph::EdgeMap<bool>& used) {
	using lemon::FullGraph;
	FullGraph::Node previous = lemon::INVALID;
	FullGraph::Node currentCity = g(0);
	size_t indexInTour = 0;
	do {
		if (indexInTour >= order.size()) {
			throw std::runtime_error("Invalid TSP solution, includes a cycle not containing node 1");
		}
		order[indexInTour] = FullGraph::id(currentCity);
		++indexInTour;
		bool foundNext = false;
		for (FullGraph::OutArcIt it(g, currentCity); it != lemon::INVALID; ++it) {
			FullGraph::Node next = g.target(it);
			if (next != previous && next != currentCity && used[it]) {
				cost += inst->getDistance(FullGraph::id(currentCity), FullGraph::id(next));
				previous = currentCity;
				currentCity = next;
				foundNext = true;
				break;
			}
		}
		if (!foundNext) {
			throw std::runtime_error("Invalid TSP solution, includes node with degree 1");
		}
	} while (currentCity != g(0));
	if (indexInTour < order.size()) {
		throw std::runtime_error("Invalid TSP solution, node 1 is in a short cycle");
	}
}

TSPSolution::TSPSolution(const TSPInstance& inst, const lemon::FullGraph& g,
						 const lemon::FullGraph::EdgeMap<bool>& used)
		: inst(&inst), order(static_cast<size_t>(inst.getCityCount())) {
	initFromGraph(g, used);
}

TSPSolution::TSPSolution(const TSPInstance& inst, std::vector<city_id> order)
		: inst(&inst), order(std::move(order)) {
	initTourCost();
}

/**
 * Liest eine Tour im TSPLib-Format ein
 * @param instance Die zur Tour gehörende Instanz
 */
TSPSolution::TSPSolution(const TSPInstance& instance, std::istream& input) : inst(&instance) {
	std::string line;
	bool emptyLines = false;
	while (std::getline(input, line)) {
		if (!line.empty()) {
			if (emptyLines) {
				std::cout << "Skipped empty line(s)" << std::endl;
			}
			std::stringstream ss(line);
			std::string keyword = tsp_util::readKeyword(ss);
			if (keyword == "NAME" || keyword == "COMMENT") {
				//NOP
			} else if (keyword == "TYPE") {
				std::string type;
				ss >> type;
				if (type != "TOUR") {
					throw std::runtime_error("Tried to read tour from non-tour file");
				}
			} else if (keyword == "DIMENSION") {
				city_id dim;
				ss >> dim;
				if (dim != inst->getCityCount()) {
					throw std::runtime_error("Tour and instance have different dimension");
				}
			} else if (keyword == "TOUR_SECTION") {
				order.resize(inst->getCityCount());
				for (city_id& i : order) {
					i = tsp_util::readOrThrow<city_id>(input) - 1;
					if (i < 0) {
						throw std::runtime_error("Invalid city in tour: "
												 + std::to_string(i));
					}
				}
				auto end = tsp_util::readOrThrow<city_id>(input);
				if (end != -1) {
					throw std::runtime_error("Tour wasn't followed by -1!");
				}
			} else if (keyword == "EOF") {
				break;
			} else {
				throw std::runtime_error("Unknown keyword in tour file: " + keyword);
			}
		} else {
			emptyLines = true;
		}
	}
	if (order.empty()) {
		throw std::runtime_error("Tour file did not contain tour data!");
	}
	initTourCost();
}

/**
 * @return Die Kosten der Lösung
 */
cost_t TSPSolution::getCost() const {
	return cost;
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
	if (inst == nullptr) {
		throw std::runtime_error("Tried to write invalid TSP solution!");
	}
	out << "NAME: " << inst->getName() << ".tsp.tour\n"
		<< "TYPE: TOUR\n"
		<< "DIMENSION: " << inst->getCityCount() << "\n"
		<< "TOUR_SECTION\n";
	for (city_id c:getOrder()) {
		out << c + 1 << "\n";
	}
	out << "-1\n";
}

bool TSPSolution::isValid() const {
	return inst != nullptr;
}

/**
 * Berechnet die Kosten der Tour aus der Knotenreihenfolge
 */
void TSPSolution::initTourCost() {
	cost = 0;
	for (city_id i = 0; i < inst->getCityCount(); ++i) {
		city_id cityA = order[i];
		city_id cityB = order[(i + 1) % inst->getCityCount()];
		cost += inst->getDistance(cityA, cityB);
	}
}
