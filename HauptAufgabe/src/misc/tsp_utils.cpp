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
			if (tolerance.positive(solution[varId])) {
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
