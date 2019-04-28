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
			/*{"burma14",   3323},
			{"ulysses16", 6859},
			{"gr17",      2085},
			{"gr21",      2707},
			{"ulysses22", 7013},
			{"gr24",      1272},
			{"fri26",     937},
			{"bays29",    2020},
			{"bayg29",    1610},
			{"swiss42",   1273},
			{"dantzig42", 699},
			{"gr48",      5046},
			{"att48",     10628},
			{"hk48",      11461},
			{"eil51",     426},
			{"berlin52",  7542},
			{"brazil58",  25395},
			{"st70",      675},
			{"eil76",     538},
			{"pr76",      108159},
			{"gr96",      55209},
			{"rat99",     1211},
			{"kroA100",   21282},
			{"kroB100",   22141},
			{"kroC100",   20749},
			{"kroD100",   21294},
			{"kroE100",   22068},
			{"rd100",     7910},
			{"eil101",    629},
			{"lin105",    14379},
			{"pr107",     44303},
			{"gr120",     6942},
			{"pr124",     59030},
			{"bier127",   118282},
			{"ch130",     6110},
			{"pr136",     96772},
			{"gr137",     69853},
			{"pr144",     58537},
			{"kroA150",   26524},
			{"kroB150",   26130},
			{"ch150",     6528},
			{"pr152",     73682},
			{"u159",      42080},
			{"si175",     21407},
			{"brg180",    1950},
			{"rat195",    2323},
			{"d198",      15780},
			{"kroA200", 29368},
			{"kroB200", 29437},
			{"gr202",   40160},
			//{"ts225",     126643},//Braucht aus mir nicht ersichtlichen Gründen sehr lange (>10h)
			{"tsp225",  3916},
			{"pr226",   80369},
			{"gr229",   134602},
			{"gil262",  2378},
			{"pr264",   49135},
			{"a280",    2579},
			{"pr299",   48191},
			//{"lin318",    42029},//Hat eine "festgelegte" Kante, weil es ursprünglich ein Hamilton-Pfad war
			//{"linhp318",  41345},//Siehe ^
			{"rd400",   15281},
			{"fl417",   11861},
			{"gr431",   171414},
			{"pr439",   107217},
			{"pcb442",  50778},
			{"d493",    35002},
			{"att532",  27686},
			{"si535",   48450},
			{"ali535",  202339},
			{"pa561",   2763},
			{"u574",    36905},
			{"rat575",  6773},
			{"p654",    34643},
			{"d657",    48913},//Ist auf der Website falsch
			{"gr666",   294358},*/
			{"u724",    41910},
			{"rat783",  8806},
			{"dsj1000", 18660188},
			{"pr1002",  259045},
			{"si1032",  92650},
			{"u1060",   224094},
			{"vm1084",  239297},
			{"pcb1173", 56892},
			{"d1291",   50801},
			{"rl1304",  252948},
			{"rl1323",  270199},
			{"nrw1379", 56638},
			//{"fl1400",  20127},
			{"u1432",   152970},
			//{"fl1577",  22249},
			{"d1655",   62128},
			{"vm1748",  336556},
			{"u1817",   57201},
			{"rl1889",  316536},
			{"d2103",   80450},
			{"u2152",   64253},
			{"u2319",   234256},
			{"pr2392",  378032},
	};
	int status;
	CPXENVptr env = CPXopenCPLEX(&status);
	if (status != 0) {
		throw std::runtime_error("Failed to open CPLEX environment: " + std::to_string(status));
	}
	std::vector<double> times;
	for (const auto& inst:instances) {
		clock_t start = std::clock();
		std::cout << "Solving instance " << inst.first << std::endl;
		std::ifstream in("../instances/" + inst.first + ".tsp");
		if (!in) {
			std::cout << "Instance not found!" << std::endl;
			continue;
		}
		TSPInstance tsp(in);
		TSPSolution init = tspsolvers::solveGreedy(tsp);
		std::cout << "Using initial solution with cost " << init.getCost() << std::endl;
		TSPSolution sol = tspsolvers::solveLP(tsp, &init, env, false);
		clock_t end = std::clock();
		double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
		if (sol.getCost() != inst.second) {
			std::cerr << "Failed to solve " << inst.first << ", found solution has cost " << sol.getCost()
					  << ", but optimal cost is " << inst.second << std::endl;
			break;
		} else {
			std::cout << "Found correct solution for " << inst.first << " in " << elapsed_secs << " seconds"
					  << std::endl;
		}
		times.push_back(elapsed_secs);
	}
	if (instances.size() == times.size()) {
		for (size_t i = 0; i < instances.size(); ++i) {
			std::cout << instances[i].first << ": " << times[i] << std::endl;
		}
	}
	CPXcloseCPLEX(&env);
}
