#include <string>
#include <fstream>
#include <tsp_instance.hpp>
#include <subtour_cut_gen.hpp>
#include <two_matching_cut_gen.hpp>
#include <branch_and_cut.hpp>
#include <ctime>
#include <tsp_solution.hpp>
#include <tsp_solvers.hpp>

int main() {
	const std::vector<std::pair<std::string, cost_t>> instances{
			{"berlin52",  7542},
			//{"KorteVygenNonInt", 10},
			{"fri26",     937},
			{"dantzig42", 699},
			{"bayg29",    1610},
			{"bays29",    2020},
			{"brazil58",  25395},
			{"hk48",      11461},
			{"lin105",    14379},
			{"rd100",     7910},
			{"kroA100",   21282},
			{"kroB100",   22141},
			{"kroC100",   20749},
			{"kroD100",   21294},
			{"kroE100",   22068},
			{"kroA150",   26524},
			{"a280",      2579},
			{"bier127",   118282},
			{"eil51",     426},
			{"eil76",     538},
			{"eil101",    629},
			{"ch130",     6110},
			{"rat99",     1211},
			{"gr17",      2085},
			{"gr21",      2707},
			{"gr24",      1272},
			{"gr48",      5046},
			{"gr96",      55209},
			{"gr120",     6942},
			{"gr137",     69853},
			{"gr202",     40160},
			{"gr229",     134602},
			{"kroA200",   29368},
			{"pr299",     48191},
			{"gil262",    2378},
			{"lin318",    42029},
			{"Tnm52",     551609},
			{"gr431",     171414},
			{"fl417",     11861},
			{"pcb442",    50778},
	};
	int status;
	CPXENVptr env = CPXopenCPLEX(&status);
	if (status!=0) {
		throw std::runtime_error("Failed to open CPLEX environment: "+std::to_string(status));
	}
	for (const auto& inst:instances) {
		clock_t start = std::clock();
		std::cout << "Solving instance " << inst.first << std::endl;
		std::ifstream in("../instances/"+inst.first+".tsp");
		if (!in) {
			std::cout << "Instance not found!" << std::endl;
			continue;
		}
		TSPInstance tsp(in);
		TSPSolution init = tspsolvers::solveGreedy(tsp);
		std::cout << "Using initial solution with cost " << init.getCost() << std::endl;
		TSPSolution sol = tspsolvers::solveLP(tsp, &init, env, 1536*1024*1024);
		clock_t end = std::clock();
		double elapsed_secs = double(end-start)/CLOCKS_PER_SEC;
		if (sol.getCost()!=inst.second) {
			std::cerr << "Failed to solve " << inst.first << ", found solution has cost " << sol.getCost()
					  << ", but optimal cost is " << inst.second << std::endl;
			break;
		} else {
			std::cout << "Found correct solution for " << inst.first << " in " << elapsed_secs << " seconds"
					  << std::endl;
		}
	}
	CPXcloseCPLEX(&env);
}
