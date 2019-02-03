
#include <string>
#include <fstream>
#include <tsp_instance.hpp>
#include <subtour_cut_gen.hpp>
#include <branch_and_cut.hpp>
#include <ctime>

int main() {
	const std::vector<std::pair<std::string, cost_t>> instances{
			{"berlin52",         7542},
			{"KorteVygenNonInt", 10},
			{"fri26",            937},
			{"dantzig42",        699},
			{"bayg29",           1610},
			{"bays29",           2020},
			{"brazil58",         25395},
			{"hk48",             11461},
			{"lin105",           14379},
			{"rd100",            7910},
			{"kroA100",          21282},
			{"kroB100",          22141},
			{"kroC100",          20749},
			{"kroD100",          21294},
			{"kroE100",          22068},
	};
	for (const auto& inst:instances) {
		clock_t start = std::clock();
		std::cout << "Solving instance " << inst.first << std::endl;
		std::ifstream in("../instances/"+inst.first+".tsp");
		if (!in) {
			std::cout << "Instance not found!" << std::endl;
			continue;
		}
		TSPInstance tsp(in);
		LinearProgram lp(inst.first, LinearProgram::minimize);
		tsp.setupBasicLP(lp);
		SubtourCutGen subtours(tsp);
		BranchAndCut bac(lp, {&subtours});
		std::vector<long> tour = bac.solve();
		cost_t totalCost = 0;
		for (variable_id varId = 0; varId<tour.size(); ++varId) {
			if (tour[varId]>0) {
				Graph::Edge e = tsp.getEdge(varId);
				totalCost += tsp.getGraphDistances()[e];
			}
		}
		clock_t end = std::clock();
		double elapsed_secs = double(end-start)/CLOCKS_PER_SEC;
		std::cout << "Found solution of length " << totalCost << " (Optimum: " << inst.second << ")" << " for "
				  << inst.first
				  << " in " << elapsed_secs << " seconds"
				  << std::endl;
	}
}