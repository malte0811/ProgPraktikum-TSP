#include <utility>
#include "tour.hpp"

Tour::Tour(std::vector<Vec2> nodes, std::vector<node_id> order, EdgeWeightType eT,
		   cost_t length): nodeLocations(std::move(nodes)), order(std::move(order)), edgeType(eT), length(length) {}
