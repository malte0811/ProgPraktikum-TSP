#include <two_matching_cut_gen.hpp>
#include <lemon/unionfind.h>
#include <union_find.hpp>
#include <tsp_utils.hpp>
#include <tsp_lp_data.hpp>
#include <blossom_finder.hpp>

TwoMatchingCutGen::TwoMatchingCutGen(const TSPInstance& inst, const TspLpData& lpData, bool contract)
		: tsp(inst), enableContraction(contract), tolerance(1e-5), lpData(lpData) {}

/**
 * Findet verletzte 2-Matching-Constraints, falls es solche gibt. Der Algorithmus entspricht grob dem aus Satz 12.21 mit
 * einigen Änderungen:
 * 1. z ist hier nicht nötig: Die Grad-Constraints sind immer erfüllt, also gilt für alle an z anliegenden Kanten e
 * c(e)=0 und c'(e)=inf. Sie können also nie in der Menge F einer verletzten 2-Matching-Constraint enthalten sein.
 * 2. Der Graph wird vor dem Anwenden von Lemma 12.20 durch Pfad-Kontraktion nach Proposition 4.6 in
 * Padberg, M. & Rinaldi, G. Mathematical Programming (1990) 47: 219. https://doi.org/10.1007/BF01580861
 * 3. Im Algorithmus aus Lemma 12.21 werden hier nicht die minimalen Blüten berechnet, sondern alle aus dem Gomory-Hu-
 * Baum entstehenden Blüten mit Wert echt kleiner 1.
 * 4. Die Blüten werden vor dem Hinzufügen der Constraints so verändert, dass F ein Matching ist. Der Algorithmus wird
 * im selben Paper wie oben ohne Beweis gegeben. TODO: Quelle mit Beweis finden oder selbst einen schreiben
 */
CutGenerator::CutStatus TwoMatchingCutGen::validate(LinearProgram& lp, const std::vector<double>& solution,
													CutStatus currentStatus) {
	if (currentStatus == CutGenerator::maybe_recalc) {
		return CutGenerator::valid;
	}
	Graph workGraph;
	//Ordnet den Städten der TSP-Instanz einen Knoten im Arbeitsgraphen zu. Nur bis zum Aufruf von contractPaths gültig!
	std::vector<Graph::Node> origToWork(static_cast<size_t>(tsp.getCityCount()));
	//c wie in Satz 12.21
	Graph::EdgeMap<double> c(workGraph);
	Graph::NodeMap <city_id> workToOrig(workGraph);
	tsp_util::addSupportGraphEdges(tsp, lpData, tolerance, solution, workGraph,
								   origToWork, workToOrig, c);
	BlossomFinder finder(workGraph, c, tolerance, enableContraction);
	std::vector<BlossomFinder::Blossom> allMin = finder.findViolatedBlossoms();
	if (!allMin.empty()) {
		std::set<tsp_util::ConstraintWithSlack, tsp_util::CompareOrderedConstraint> allConstrs;
		for (BlossomFinder::Blossom& min:allMin) {
			std::vector<bool> inHandle(tsp.getCityCount());
			for (Graph::Node n:min.handle) {
				inHandle[workToOrig[n]] = true;
			}
			for (Graph::Edge e:min.teeth) {
				Graph::Node u = workGraph.u(e);
				Graph::Node v = workGraph.v(e);
				assert(inHandle[workToOrig[u]] != inHandle[workToOrig[v]]);
			}

			const size_t sizeF = min.teeth.size();
			/*
			 * Mit |F|==1 ist auch die Subtour-Constraint für X verletzt und impliziert die
			 * 2-Matching-Constraint
			 */
			assert(sizeF%2==1);
			//Die Indizes der Variablen in den hinzugefügten Constraints
			std::vector<variable_id> indices;
			std::vector<double> coeffs(lpData.getVariableCount());
			if (sizeF>1) {
				for (Graph::Edge e:min.teeth) {
					Graph::Node u = workGraph.u(e);
					Graph::Node v = workGraph.v(e);
					variable_id var = lpData.getVariable(workToOrig[u], workToOrig[v]);
					TspLpData::Edge e1 = lpData.getEdge(var);
					assert(inHandle[e1.first] != inHandle[e1.second]);
					indices.push_back(var);
					coeffs[var] += 1;
				}
			}
			std::vector<city_id> handleTSP;
			handleTSP.reserve(min.handle.size());
			for (Graph::Node n:min.handle) {
				handleTSP.push_back(workToOrig[n]);
			}
			lpData.sparserInducedSum(handleTSP, coeffs, indices);
			if (!indices.empty()) {
				LinearProgram::Constraint constr;
				std::sort(indices.begin(), indices.end());
				if (sizeF > 1) {
					//2-Matching-Constraint
					constr = LinearProgram::Constraint::fromDense(indices, coeffs,
													   LinearProgram::less_eq,
													   static_cast<size_t>(min.handle.size() + sizeF / 2));
				} else {
					//Subtour-Constraint
					constr = LinearProgram::Constraint::fromDense(indices, coeffs,
													   LinearProgram::less_eq,
										 min.handle.size() - 1);
				}
				double lhs = constr.evalLHS(solution);
				if (tolerance.less(constr.getRHS(), lhs)) {
					allConstrs.insert({constr, constr.getRHS()-lhs});
				}
			}
		}
		if (!allConstrs.empty()) {
			size_t maxNonzero = 0, sumNZ = 0;
			std::vector<LinearProgram::Constraint> toAdd;
			for (const tsp_util::ConstraintWithSlack& cws:allConstrs) {
				toAdd.push_back(cws.constraint);
				const size_t nz = cws.constraint.getNonzeroes().size();
				sumNZ += nz;
				if (nz>maxNonzero) {
					maxNonzero = nz;
				}
				if (sumNZ>(tsp.getCityCount()*tsp.getCityCount())) {
					break;
				}
			}
			std::cout << "2Matching: Max nonzero: " << maxNonzero << ", sum: " << sumNZ << std::endl;
			lp.addConstraints(toAdd);
			return CutGenerator::maybe_recalc;
		}
	}
	return CutGenerator::valid;
}