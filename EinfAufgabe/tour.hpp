#ifndef TOUR_HPP
#define TOUR_HPP

#include <vector>
#include <array>
#include <ostream>
#include "graph.hpp"

class Tour {
public:
	using node_id = Graph::node_id;
	using edge_id = Graph::edge_id;
	using cost_t = Graph::cost_t;
	using Vec2 = Graph::Vec2;

	Tour(std::vector<node_id> order, cost_t length, std::string name);

	const cost_t length;
	const std::string name;
	const std::vector<node_id> order;

	void print(std::ostream& out);
};


#endif //TOUR_HPP
