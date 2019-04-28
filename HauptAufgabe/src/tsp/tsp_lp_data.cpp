#include <tsp_lp_data.hpp>
#include <cstddef>
#include <limits>
#include <memory>
#include <tsp_utils.hpp>
#include <tsp_instance.hpp>
#include <tsp_solution.hpp>
#include <branch_and_cut.hpp>
#include <linear_program.hpp>

using std::size_t;

TspLpData::TspLpData(const TSPInstance& inst, const TSPSolution *initial)
		: inst(inst), variableToEdge(inst.getEdgeCount()), edgeToVariable(inst.getCityCount() - 1),
		  removalBound(inst.getEdgeCount(), -std::numeric_limits<double>::max()), tolerance(0.1) {
	//Allen Kanten im vollständigen Graphen Variablen zuordnen
	variable_id currVar = 0;
	for (city_id higher = 1; higher < inst.getCityCount(); ++higher) {
		edgeToVariable[higher - 1].resize(higher);
		for (city_id lower = 0; lower < higher; ++lower) {
			variableToEdge[currVar] = {lower, higher};
			edgeToVariable[higher - 1][lower] = currVar;
			++currVar;
		}
	}
	if (initial != nullptr) {
		upperBound = *initial;
	}
}

/**
 * Setzt die beste bekannte Schranke auf die variables entsprechende Tour und entfernt dadurch unnötig gewordene
 * Variablen
 * @param variables Die LP-Variablen für die neue obere Schranke
 * @return die entfernten Variablen
 */
std::vector<variable_id> TspLpData::removeOnUpperBound(const std::vector<value_t>& variables) {
	std::vector<bool> asBools(variables.size());
	cost_t bound = 0;
	for (size_t i = 0; i < asBools.size(); ++i) {
		if (variables[i] != 0) {
			asBools[i] = true;
			bound += getCost(i);
		}
	}
	if (upperBound.isValid() && bound >= upperBound.getCost()) {
		return {};
	}
	upperBound = TSPSolution(inst, asBools, *this);
	return removeVariables();

}

/**
 * Verbessert die "Entfernungsschranken" auf die aus der LP-Lösung und den reduzierten Kosten erhaltenen Werte und
 * entfernt alle dadurch unnötig gewordenen Variablen
 * @param rootSol Eine Lösung des Wurzelknotens im B&B-Suchbaum (die nicht unbedingt alle Cut-Generatoren erfüllt)
 * @return die entfernten Variablen
 */
std::vector<variable_id> TspLpData::removeOnRootSolution(const LinearProgram::Solution& rootSol) {
	for (variable_id i = 0; i < getVariableCount(); ++i) {
		/*
		 * Entfernen, falls L+red>U-1 <=> U<L+1+red (U: obere Schranke, L: untere Schranke aus dem LP, red: reduzierte
		 * Kosten)
		 */
		double lpThreshold = rootSol.getValue() + 1 + rootSol.getReducedCosts()[i];
		if (lpThreshold > removalBound[i] && rootSol[i] == 0) {
			removalBound[i] = lpThreshold;
		}
	}
	if (!upperBound.isValid()) {
		//Es gibt noch keine obere Schranke
		return {};
	} else {
		return removeVariables();
	}
}

/**
 * Entfernt alle Variablen, die in einer Lösung, die um min. 1 besser als die beste bekannte, sicher Wert 0 haben.
 * @return die entfernten Variablen
 */
std::vector<variable_id> TspLpData::removeVariables() {
	const cost_t bound = upperBound.getCost();
	std::vector<variable_id> toRemove;
	//Die nächste freie neue Variablen-ID
	variable_id newId = 0;
	for (variable_id oldId = 0; oldId < getVariableCount(); ++oldId) {
		Edge e = variableToEdge[oldId];
		if (tolerance.less(bound, removalBound[oldId])) {
			toRemove.push_back(oldId);
			edgeToVariable[e.second - 1][e.first] = LinearProgram::invalid_variable;
		} else {
			edgeToVariable[e.second - 1][e.first] = newId;
			++newId;
		}
	}
	tsp_util::eraseEntries(removalBound, toRemove);
	tsp_util::eraseEntries(variableToEdge, toRemove);
	return toRemove;
}

/**
 * Fügt Variablen für alle Kanten und Gradbedingungen für alle Knoten zum LP hinzu
 * @param lp das zu initialisierende LP
 */
void TspLpData::setupBasicLP(LinearProgram& lp) const {
	variable_id varCount = getVariableCount();
	std::vector<double> objCoeffs(varCount);
	for (variable_id i = 0; i < varCount; ++i) {
		objCoeffs[i] = getCost(i);
	}
	std::vector<double> lower(varCount, 0);
	std::vector<double> upper(varCount, 1);
	lp.addVariables(objCoeffs, lower, upper);
	std::vector<LinearProgram::Constraint> constrs;
	constrs.reserve(inst.getCityCount());
	for (city_id i = 0; i < inst.getCityCount(); ++i) {
		std::vector<variable_id> indices;
		for (city_id otherEnd = 0; otherEnd < inst.getCityCount(); ++otherEnd) {
			if (otherEnd != i && getVariable(i, otherEnd) != LinearProgram::invalid_variable) {
				indices.push_back(getVariable(i, otherEnd));
			}
		}
		constrs.emplace_back(indices, std::vector<double>(indices.size(), 1), LinearProgram::equal,
							 2);
	}
	lp.addConstraints(constrs);
}

const TSPSolution& TspLpData::getUpperBound() const {
	return upperBound;
}

/**
 * Erhöht die Koeffizienten der induzierten Kanten um 1.
 * @param set Die induzierende Menge
 * @param values Die (dichten) Koeffizienten der Constraint
 * @param usedVars Die Variablen mit Wert ungleich 0
 * @return Die RHS der zur verwendeten inuzierten Menge gehörenden Subtour-Constraint (set.size()-1)
 */
city_id TspLpData::inducedSum(const std::vector<city_id>& set, std::vector<double>& values,
							  std::vector<variable_id>& usedVars) const {
	for (size_t i = 1; i < set.size(); ++i) {
		for (size_t j = 0; j < i; ++j) {
			variable_id var = getVariable(set[i], set[j]);
			if (var != LinearProgram::invalid_variable) {
				if (values[var] == 0) {
					usedVars.push_back(var);
				}
				values[var] += 1;
			}
		}
	}
	return set.size() - 1;
}

/**
 * Erhöht die Koeffizienten der induzierten Kanten um 1. Falls dadurch eine dünner besetzte Constraint entsteht, werden
 * die vom Komplement induzierten Kanten verwendet (x(X)<=|X|-k <=> x(V\X)<=|V\X|-k).
 * @param set Die induzierende Menge
 * @param values Die (dichten) Koeffizienten der Constraint
 * @param usedVars Die Variablen mit Wert ungleich 0
 * @return Die RHS der zur verwendeten inuzierten Menge gehörenden Subtour-Constraint
 */
city_id TspLpData::sparserInducedSum(const std::vector<city_id>& set, std::vector<double>& values,
									 std::vector<variable_id>& usedVars) const {
	/*
	 * Die Menge mit weniger induzierten Kanten ist in vollstd. Graphen die mit weniger Knoten. Wenn viele Variablen
	 * entfernt wurden, kann das prinzipiell die falsche/dichtere Menge sein. Praktisch tritt dies aber so selten auf
	 * (<<1%), dass eine genauere Entscheidung mehr Zeit verbraucht als durch die dünner besetzten Constraints gespart
	 * wird
	 */
	bool useSet = static_cast<city_id>(set.size()) <= inst.getCityCount() / 2;
	if (useSet) {
		return inducedSum(set, values, usedVars);
	} else {
		std::vector<bool> inSet(inst.getCityCount());
		for (city_id c:set) {
			inSet[c] = true;
		}
		std::vector<city_id> complement;
		complement.reserve(inst.getCityCount() - set.size());
		for (city_id i = 0; i < inst.getCityCount(); ++i) {
			if (!inSet[i]) {
				complement.push_back(i);
			}
		}
		return inducedSum(complement, values, usedVars);
	}
}
