#ifndef TOUR_HPP
#define TOUR_HPP
#include <vector>
#include <array>

enum EdgeWeightType {
	euc_2d,
	ceil_2d
};
class Tour {
public:
	//Aus graph.hpp kopiert. TODO schönere Alternative
	using node_id = unsigned;
	using edge_id = unsigned;
	using cost_t = unsigned;
	using Vec2 = std::array<double, 2>;
	Tour(std::vector<Vec2> nodes, std::vector<node_id> order,
			EdgeWeightType eT, cost_t length);

	const cost_t length;
	const std::vector<Vec2> nodeLocations;
	const std::vector<node_id> order;
	const EdgeWeightType edgeType;
};


#endif //TOUR_HPP
