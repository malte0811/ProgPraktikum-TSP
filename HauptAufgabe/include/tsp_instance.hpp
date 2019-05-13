#ifndef TSP_INSTANCE_HPP
#define TSP_INSTANCE_HPP

#include <lemon/list_graph.h>
#include <istream>
#include <linear_program.hpp>
#include <string>
#include <vector>

using cost_t = unsigned long;
//Kanten müssen als int's gespeichert werden können (wegen CPLEX), also müssen die Städte erst recht in int passen
using city_id = int;
using Graph = lemon::ListGraph;

class TSPInstance {
public:
	explicit TSPInstance(std::istream& in);

	inline cost_t getDistance(city_id a, city_id b) const;

	inline city_id getCityCount() const;

	inline variable_id getEdgeCount() const;

	std::string getName() const;

	static const city_id invalid_city;

private:
	enum EdgeWeightType {
		euc_2d,
		ceil_2d,
		explicit_,
		geo,
		att
	};
	enum EdgeFormat {
		function,
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
	return static_cast<city_id>(distances.size() + 1);
}

variable_id TSPInstance::getEdgeCount() const {
	city_id size = getCityCount();
	return (size * (size - 1)) / 2;
}

cost_t TSPInstance::getDistance(city_id a, city_id b) const {
   //TODO: Hier ist es generell eine gute Idee einen rangecheck einzubauen, ggf. nur im DEBUG modus.
	if (a > b) {
		return distances[a - 1][b];
	} else if (b > a) {
		return distances[b - 1][a];
	} else {
		return 0;
	}
}

#endif