#ifndef TSP_UTILS_HPP
#define TSP_UTILS_HPP

#include <lemon/tolerance.h>
#include <iostream>
#include <linear_program.hpp>
#include <stdexcept>
#include <string>
#include <tsp_instance.hpp>
#include <vector>

class TspLpData;

namespace tsp_util {
	using ContractionMapTSP = Graph::NodeMap<std::vector<city_id>>;
	using ContractionMapGraph = Graph::NodeMap<std::vector<Graph::Node>>;

	struct ConstraintWithSlack {
		ConstraintWithSlack() = delete;

		explicit operator const LinearProgram::Constraint&() const;

		explicit operator LinearProgram::Constraint&();

		bool operator<(const ConstraintWithSlack& other) const;

		LinearProgram::Constraint constraint;
		double slack;
	};

	void addSupportGraphEdges(const TspLpData& lpData, lemon::Tolerance<double> tolerance,
							  const std::vector<double>& solution, Graph& workGraph,
							  std::vector<Graph::Node>& origToWork, Graph::NodeMap <city_id>& workToOrig,
							  Graph::EdgeMap<double>& c);

	template<typename T>
	void eraseEntries(std::vector<T>& vec, const std::vector<variable_id>& toRemove);

	std::string readKeyword(std::istream& in);

	/**
	 * Liest einen Wert aus dem gegebenen istrean und wirft einen Fehler, falls dies nicht möglich ist
	 * @tparam T Der Typ des zu lesenden Wertes
	 */
	template<typename T>
	T readOrThrow(std::istream& input);
}


template<typename T>
void tsp_util::eraseEntries(std::vector<T>& vec, const std::vector<variable_id>& toRemove) {
	if (toRemove.empty()) {
		return;
	}
	std::vector<T> result;
	result.reserve(vec.size() - toRemove.size());
	variable_id removeIndex = 0;
	for (variable_id i = 0; i < static_cast<variable_id>(vec.size()); ++i) {
		if (removeIndex < static_cast<variable_id>(toRemove.size()) && toRemove[removeIndex] == i) {
			++removeIndex;
		} else {
			result.push_back(vec[i]);
		}
	}
	vec = result;
}

template<typename T>
T tsp_util::readOrThrow(std::istream& input) {
	T ret;
	input >> ret;
	if (!input) {
		throw std::runtime_error("Could not read input!");
	}
	return ret;
}


#endif