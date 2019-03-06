#ifndef TSP_INSTANCE_HPP
#define TSP_INSTANCE_HPP

#include <istream>
#include <vector>
#include <lemon/list_graph.h>
#include <linear_program.hpp>
#include <cmath>

using cost_t = unsigned;
using city_id = int;
using Graph = lemon::ListGraph;

class TSPInstance {
public:
	explicit TSPInstance(std::istream& in);

	inline cost_t getDistance(city_id a, city_id b) const;

	inline cost_t getCost(variable_id e) const;

	inline city_id getCityCount() const;

	inline variable_id getVariable(city_id a, city_id b) const;

	inline city_id getHigherEnd(variable_id var) const;

	inline city_id getLowerEnd(variable_id var) const;

	inline variable_id getEdgeCount() const;

	void setupBasicLP(LinearProgram& lp) const;

	std::string getName() const;

	static constexpr city_id invalid_city = -1;

private:
	enum EdgeWeightType {
		euc_2d,
		ceil_2d,
		explicit_,
		geo,
		att
	};
	enum EdgeFormat {
		full_matrix,
		lower_diag_row,
		upper_diag_row,
		upper_row
	};

	void readNodes(std::istream& input, EdgeWeightType type);

	void readEdges(std::istream& input, EdgeFormat type);

	void setDistance(city_id a, city_id b, cost_t dist);

	//distances[a][b] mit a>=b gibt die Distanz zwischen den Städten a+1 und b an
	std::vector<std::vector<cost_t>> distances;
	//Der Name der Instanz
	std::string name;
};


city_id TSPInstance::getCityCount() const {
	return static_cast<city_id>(distances.size()+1);
}

variable_id TSPInstance::getVariable(city_id a, city_id b) const {
	if (b>a) {
		std::swap(a, b);
	}
	return (a*(a-1))/2+b;
}

city_id TSPInstance::getHigherEnd(variable_id var) const {
	//TODO ist nicht wirklich schön...
	return static_cast<city_id>(.5+std::sqrt(.25+2*var));
}

city_id TSPInstance::getLowerEnd(variable_id var) const {
	city_id higher = getHigherEnd(var);
	return var-(higher*(higher-1))/2;
}

variable_id TSPInstance::getEdgeCount() const {
	city_id size = getCityCount();
	return (size*(size-1))/2;
}

cost_t TSPInstance::getDistance(city_id a, city_id b) const {
	if (a>b) {
		return distances[a-1][b];
	} else if (b>a) {
		return distances[b-1][a];
	} else {
		return 0;
	}
}

cost_t TSPInstance::getCost(variable_id e) const {
	return getDistance(getHigherEnd(e), getLowerEnd(e));
}

#endif