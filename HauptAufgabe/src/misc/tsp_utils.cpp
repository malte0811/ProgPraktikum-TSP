#include <tsp_utils.hpp>
#include <tsp_lp_data.hpp>

std::vector<variable_id> tsp_util::createFractionalGraph(const TSPInstance& tsp, const TspLpData& lpData,
														 lemon::Tolerance<double> tolerance,
														 const std::vector<double>& solution, Graph& workGraph,
														 Graph::NodeMap <city_id>& workToOrig,
														 std::vector<Graph::Node>& origToWork,
														 Graph::EdgeMap <variable_id>& toVariable,
														 Graph::EdgeMap<double>& c, Graph::NodeMap<bool>& odd) {
	for (city_id i = 0; i<tsp.getCityCount(); ++i) {
		Graph::Node newNode = workGraph.addNode();
		origToWork[i] = newNode;
		workToOrig[newNode] = {i};
	}
	std::vector<variable_id> oneEdges;
	for (city_id lower = 0; lower < tsp.getCityCount() - 1; ++lower) {
		for (city_id higher = lower + 1; higher < tsp.getCityCount(); ++higher) {
			variable_id varId = lpData.getVariable(higher, lower);
			/*
			 * Kanten mit Wert 0 können nicht in F enthalten sein, da der Wert der Blüte dann schon mindestens 1 wäre
			 * (und damit die Constraint nicht verletzt ist)
			 */
			if (solution[varId] > 0) {
				Graph::Node endU = origToWork[lower];
				Graph::Node endV = origToWork[higher];
				if (tolerance.less(solution[varId], 1)) {
					//Kante hinzufügen
					Graph::Edge eWork = workGraph.addEdge(endU, endV);
					toVariable[eWork] = varId;
					c[eWork] = solution[varId];
				} else {
					/*
					 * Kanten im Schnitt mit Wert 1 müssen in F enthalten sein, sonst ist der Wert mindestens 1. Durch
					 * das Hinzufügen bzw Entfernen der Enden aus T bleiben die möglichen Blüten gleich.
					 */
					odd[endU] = !odd[endU];
					odd[endV] = !odd[endV];
					oneEdges.push_back(varId);
				}
			}
		}
	}
	return oneEdges;
}

void tsp_util::addSupportGraphEdges(const TSPInstance& tsp, const TspLpData& lpData, lemon::Tolerance<double> tolerance,
									const std::vector<double>& solution, Graph& workGraph,
									std::vector<Graph::Node>& origToWork,
									Graph::NodeMap <city_id>& workToOrig, Graph::EdgeMap<double>& c) {
	if (Graph::NodeIt(workGraph) == lemon::INVALID) {
		origToWork.resize(tsp.getCityCount());
		for (city_id i = 0; i < tsp.getCityCount(); ++i) {
			Graph::Node newNode = workGraph.addNode();
			workToOrig[newNode] = i;
			origToWork[i] = newNode;
		}
	} else {
		//Bereits existierende Kanten entfernen
		for (Graph::EdgeIt it(workGraph); it != lemon::INVALID;) {
			Graph::Edge e = it;
			++it;
			workGraph.erase(e);
		}
	}
	/*
	 * Alle Kanten einfügen, deren Variablen einen echt positiven Wert haben
	 */
	for (variable_id i = 0; i < solution.size(); ++i) {
		if (tolerance.positive(solution[i])) {
			TspLpData::Edge e = lpData.getEdge(i);
			Graph::Edge inWork = workGraph.addEdge(origToWork[e.first], origToWork[e.second]);
			c[inWork] = solution[i];
		}
	}
}

std::string tsp_util::readKeyword(std::istream& in) {
	std::string keyword;
	in >> keyword >> std::ws;
	/*
	 * Doppelpunkt ignorieren. Der Doppelpunkt kann sowohl direkt hinter dem "Operator" als auch durch Whitespace
	 * abgetrennt auftreten.
	 */
	if (in.peek() == ':') {
		in.ignore();
	} else if (keyword.back() == ':') {
		keyword = keyword.substr(0, keyword.size() - 1);
	}
	return keyword;
}
