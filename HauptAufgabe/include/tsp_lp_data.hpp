#ifndef TSP_LP_DATA_HPP
#define TSP_LP_DATA_HPP

#include <lemon/tolerance.h>
#include <tsp_instance.hpp>
#include <tsp_solution.hpp>
#include <utility>
#include <variable_remover.hpp>
#include <vector>
#include <branch_and_cut.hpp>
#include <linear_program.hpp>

/*
 * Speichert, welche LP-Variablen welchen Kanten entsprechen, entscheidet, welche Variablen entfernt werden können und
 * stellt einige allgemeine Funktionen zum TSP-LP zur Verfügung
 */
class TspLpData : public VariableRemover {
public:
	using Edge = std::pair<city_id, city_id>;

	TspLpData(const TSPInstance& inst, const TSPSolution *initial);

	std::vector<variable_id> removeOnUpperBound(const std::vector<value_t>& variables) override;

	std::vector<variable_id> removeOnRootSolution(const LinearProgram::Solution& rootSol) override;

	std::vector<variable_id> removeVariables();

	inline cost_t getCost(variable_id var) const;

	inline Edge getEdge(variable_id var) const;

	void setupBasicLP(LinearProgram& lp) const;

	inline variable_id getVariable(city_id a, city_id b) const;

	const TSPSolution& getUpperBound() const;

	inline variable_id getVariableCount() const;

	inline const TSPInstance& getTSP() const;

	city_id inducedSum(const std::vector<city_id>& set, std::vector<double>& values,
					   std::vector<variable_id>& usedVars) const;

	city_id sparserInducedSum(const std::vector<city_id>& set, std::vector<double>& values,
							  std::vector<variable_id>& usedVars) const;

private:
	const TSPInstance& inst;
	std::vector<Edge> variableToEdge;
	//Ordnet einer Kante ihre Variable zu. Größerer Index zuerst.
	std::vector<std::vector<variable_id>> edgeToVariable;
	/*
	 * Speichert, wie klein die obere Schranke werden muss, damit die dem Index entsprechende Kante in keiner besseren
	 * Lösung vorkommen kann
	 */
	std::vector<double> removalBound;
	//Die aktuelle obere Schranke
	TSPSolution upperBound;
	lemon::Tolerance<double> tolerance;
};

/**
 * @return die zur gegebenen Variable gehörende Kante
 */
TspLpData::Edge TspLpData::getEdge(variable_id var) const {
	return variableToEdge[var];
}

/**
 * @return die Kosten der Variablen, d.h. die Länge der zugehörigen Kante
 */
cost_t TspLpData::getCost(variable_id var) const {
	Edge e = getEdge(var);
	return inst.getDistance(e.first, e.second);
}

/**
 * @return Die Variable für die Kante zwischen a und b, oder LinearProgram::invalid_variable, falls es keine gibt
 */
variable_id TspLpData::getVariable(city_id a, city_id b) const {
	if (a == b) {
		return LinearProgram::invalid_variable;
	}
	if (a < b) {
		std::swap(a, b);
	}
	return edgeToVariable[a - 1][b];
}

variable_id TspLpData::getVariableCount() const {
	return variableToEdge.size();
}

const TSPInstance& TspLpData::getTSP() const {
	return inst;
}

#endif