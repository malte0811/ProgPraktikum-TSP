#ifndef TSP_UTILS_HPP
#define TSP_UTILS_HPP

#include <linear_program.hpp>
#include <tsp_instance.hpp>
#include "tsp_lp_data.hpp"

namespace tsp_util {
	using ContractionMapTSP = Graph::NodeMap<std::vector<city_id>>;
	using ContractionMapGraph = Graph::NodeMap<std::vector<Graph::Node>>;
	std::vector<variable_id> createFractionalGraph(const TSPInstance& tsp, const TspLpData& lpData,
												   lemon::Tolerance<double> tolerance,
												   const std::vector<double>& solution, Graph& workGraph,
												   Graph::NodeMap <city_id>& workToOrig,
												   std::vector<Graph::Node>& origToWork,
												   Graph::EdgeMap <variable_id>& toVariable,
												   Graph::EdgeMap<double>& c, Graph::NodeMap<bool>& odd);

	void addSupportGraphEdges(const TSPInstance& tsp, const TspLpData& lpData, lemon::Tolerance<double> tolerance,
							  const std::vector<double>& solution, Graph& workGraph,
							  std::vector<Graph::Node>& origToWork,
							  Graph::NodeMap <city_id>& workToOrig, Graph::EdgeMap<double>& c);

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
	size_t removeIndex = 0;
	for (size_t i = 0; i < vec.size(); ++i) {
		if (removeIndex < toRemove.size() && toRemove[removeIndex] == i) {
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
		throw std::invalid_argument("Could not read input!");
	}
	return ret;
}


#endif