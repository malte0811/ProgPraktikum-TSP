
#include <string>
#include <fstream>
#include <tsp_instance.hpp>
#include <subtour_cut_gen.hpp>
#include <branch_and_cut.hpp>

int main() {
	const std::vector<std::pair<std::string, cost_t>> instances{
			{"berlin52",         7542},
			{"KorteVygenNonInt", 10},
			{"fri26",            937},
			{"eil51",            426}};
	for (const auto& inst:instances) {
		std::ifstream in("../instances/"+inst.first+".tsp");
		TSPInstance tsp(in);
		LinearProgram lp(inst.first, LinearProgram::minimize);
		tsp.setupBasicLP(lp);
		SubtourCutGen subtours(tsp.getGraph());
		BranchAndCut bac(lp, {&subtours});
		std::vector<long> tour = bac.solve();
		edge_id edgeCount = 0;
		cost_t totalCost = 0;
		for (edge_id eid = 0; eid<tour.size(); ++eid) {
			if (tour[eid]>0) {
				Graph::Edge e = Graph::edgeFromId(eid);
				totalCost += tsp.getGraphDistances()[e];
				++edgeCount;
			}
		}
		std::cout << inst.first << " has a solution with length" << totalCost << " (Optimum: " << inst.second << ")"
				  << std::endl;
	}
}