#include <two_matching_cut_gen.hpp>
#include <lemon/unionfind.h>
#include <union_find.hpp>
#include <tsp_utils.hpp>
#include <tsp_lp_data.hpp>
#include <blossom_finder.hpp>

TwoMatchingCutGen::TwoMatchingCutGen(const TSPInstance& inst, const TspLpData& lpData, bool contract)
		: tsp(inst), enableContraction(contract), tolerance(1e-5), lpData(lpData) {}

CutGenerator::CutStatus TwoMatchingCutGen::validate(LinearProgram& lp, const std::vector<double>& solution,
													CutStatus currentStatus) {
	//Nur erfüllt, falls der SimpleCombCutGen schon verletzte Constraints, meist 2-Matching-Constraints, gefunden hat
	if (currentStatus == CutGenerator::maybe_recalc) {
		return CutGenerator::valid;
	}
	Graph supportGraph;
	//Ordnet den Städten der TSP-Instanz einen Knoten im Arbeitsgraphen zu. Nur bis zum Aufruf von contractPaths gültig!
	std::vector<Graph::Node> origToWork(static_cast<size_t>(tsp.getCityCount()));
	//Die Kantenkapazitäten
	Graph::EdgeMap<double> capacity(supportGraph);
	//Ordnet jedem Knoten die entsprechende Stadt zu
	Graph::NodeMap <city_id> graphToTSP(supportGraph);
	tsp_util::addSupportGraphEdges(tsp, lpData, tolerance, solution, supportGraph,
								   origToWork, graphToTSP, capacity);
	//Verletzten Blüten finden
	BlossomFinder finder(supportGraph, capacity, tolerance, enableContraction);
	std::vector<BlossomFinder::Blossom> allMin = finder.findViolatedBlossoms();
	if (!allMin.empty()) {
		//Constraints zu den Blüten berechnen und nach "Verletztheit" sortiert speichern
		std::set<tsp_util::ConstraintWithSlack, tsp_util::CompareOrderedConstraint> allConstrs;
		for (BlossomFinder::Blossom& min:allMin) {
			const size_t sizeF = min.teeth.size();
			assert(sizeF%2==1);
			//Die Indizes der Variablen in den hinzugefügten Constraints
			std::vector<variable_id> indices;
			//Die Koeffizienten für alle Variablen
			std::vector<double> coeffs(lpData.getVariableCount());
			/*
			 * Mit |F|==1 ist auch die Subtour-Constraint für X verletzt und impliziert die
			 * 2-Matching-Constraint
			 */
			if (sizeF>1) {
				//Zähne zur Contraint hinzufügen
				for (Graph::Edge e:min.teeth) {
					Graph::Node u = supportGraph.u(e);
					Graph::Node v = supportGraph.v(e);
					variable_id var = lpData.getVariable(graphToTSP[u], graphToTSP[v]);
					TspLpData::Edge e1 = lpData.getEdge(var);
					indices.push_back(var);
					coeffs[var] += 1;
				}
			}
			//Die zum Griff gehörenden Städte
			std::vector<city_id> handleTSP;
			handleTSP.reserve(min.handle.size());
			for (Graph::Node n:min.handle) {
				handleTSP.push_back(graphToTSP[n]);
			}
			//Griff zur Constraint hinzufügen
			city_id handleSubtourRHS = lpData.sparserInducedSum(handleTSP, coeffs, indices);
			if (!indices.empty()) {
				LinearProgram::Constraint constr;
				//Indizes sortieren, damit eventuell doppelte Constraints erkannt werden können
				std::sort(indices.begin(), indices.end());
				if (sizeF > 1) {
					//2-Matching-Constraint
					constr = LinearProgram::Constraint::fromDense
							(indices, coeffs, LinearProgram::less_eq, static_cast<size_t>(handleSubtourRHS + 1 +
																					  sizeF / 2));
				} else {
					//Subtour-Constraint
					constr = LinearProgram::Constraint::fromDense
							(indices, coeffs, LinearProgram::less_eq, handleSubtourRHS);
				}
				double lhs = constr.evalLHS(solution);
				//Falls die Constraint knapp nicht oder nur schwach verletzt ist, liegt dies an Rundungsfehlern
				if (tolerance.less(constr.getRHS(), lhs)) {
					allConstrs.insert({constr, constr.getRHS()-lhs});
				} else if (constr.getRHS() - lhs > 0.5) {
					//Falls die Constraint aber relativ stark verletzt ist, gibt es vermutlich einen Fehler im Algorithmus
					std::cerr << "Not violated: " << lhs << " vs " << constr.getRHS() << std::endl;
				}
			}
		}
		if (!allConstrs.empty()) {
			size_t sumNZ = 0;
			/*
			 * Constraints hinzufügen (stark verletzte zuerst), bis viele Nonzeroes zum LP hinzugefügt wurden. Dann
			 * werden keine weiteren Constraints hinzugefügt, um den Speicherverbrauch möglichst klein zu halten
			 */
			std::vector<LinearProgram::Constraint> toAdd;
			for (const tsp_util::ConstraintWithSlack& cws:allConstrs) {
				toAdd.push_back(cws.constraint);
				const size_t nz = cws.constraint.getNonzeroes().size();
				sumNZ += nz;
				if (sumNZ > (tsp.getCityCount() * tsp.getCityCount()) / 4) {
					break;
				}
			}
			lp.addConstraints(toAdd);
			return CutGenerator::maybe_recalc;
		}
	}
	return CutGenerator::valid;
}