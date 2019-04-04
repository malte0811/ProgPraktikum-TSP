#include <linear_program.hpp>
#include <vector>
#include <string>
#include <tsp_instance.hpp>
#include <tsp_solvers.hpp>
#include <tsp_solution.hpp>
#include <ctime>
#include <iostream>
#include <fstream>
#include <lemon/full_graph.h>
#include <lemon/adaptors.h>
#include <lemon/kruskal.h>
#include <lemon/bfs.h>
#include <one_tree.hpp>
#include <lemon/connectivity.h>

//TODO Kanten an 0/1?
void eliminateEdges(const TSPInstance& inst, cost_t upper, CPXENVptr env) {
	using lemon::FullGraph;
	using Tree = lemon::FilterEdges<FullGraph>;
	OneTree one(inst);
	one.run();
	const Tree& tree = one.getTree();
	const FullGraph& g = one.getGraph();
	const FullGraph::EdgeMap<bool>& inTree = one.getInTree();
	const FullGraph::EdgeMap<double>& redCosts = one.getReducedCosts();
	std::cout << "Upper bound: " << upper << ", lower: " << one.getBound() << std::endl;

	size_t removed = 0;
	for (FullGraph::NodeIt it(g); it != lemon::INVALID; ++it) {
		lemon::Bfs<Tree> bfs(tree);
		bfs.run(it);
		assert(lemon::connected(tree));
		FullGraph::NodeIt it2 = it;
		++it2;
		for (; it2 != lemon::INVALID; ++it2) {
			FullGraph::Edge e = g.edge(it, it2);
			if (inTree[e]) {
				continue;
			}
			double maxCost = 0;
			FullGraph::Node current = it2;
			while (current != it) {
				FullGraph::Arc predArc = bfs.predArc(current);
				if (redCosts[predArc] > maxCost) {
					maxCost = redCosts[predArc];
				}
				assert(current != g.source(predArc));
				current = g.source(predArc);
			}
			//TODO tolerance
			if (redCosts[e] - maxCost > upper - one.getBound() - 1) {
				++removed;
			}
		}
	}
	std::cout << "Removed " << removed << " of " << inst.getEdgeCount() << " edges" << std::endl;
}

int main() {

	const std::vector<std::pair<std::string, cost_t>> instances{
			//{"berlin52",         7542},
			//{"KorteVygenNonInt", 10},
			//{"fri26",            937},
			//{"dantzig42",        699},
			//{"bayg29",           1610},
			//{"bays29",           2020},
			//{"brazil58",         25395},
			//{"hk48",             11461},
			//{"lin105",           14379},
			//{"rd100",            7910},
			//{"kroA100",          21282},
			//{"kroB100",          22141},
			//{"kroC100",          20749},
			//{"kroD100",          21294},
			//{"kroE100",          22068},
			//{"kroA150",          26524},
			//{"a280",             2579},
			//{"bier127",          118282},
			//{"eil51",            426},
			//{"eil76",            538},
			//{"eil101",           629},
			//{"ch130",            6110},
			//{"rat99",            1211},
			//{"gr17",             2085},
			//{"gr21",             2707},
			//{"gr24",             1272},
			//{"gr48",             5046},
			//{"gr96",             55209},
			//{"gr120",            6942},
			//{"gr137",            69853},
			//{"gr202",            40160},
			//{"gr229",            134602},
			//{"kroA200",          29368},
			//{"pr299",            48191},
			//{"gil262",           2378},
			//{"lin318",           42029},
			//{"Tnm52",            551609},
			//{"pcb442",           50778},
			{"d493", 35002},//52634 of 121278
			{"fl417", 11861},//41739 of 86736
			{"d657", 48912},//117671 of 215496
			{"gr666", 294358},//52732 of 221445
			{"nrw1379", 56638},//486849 of 950131
	};
	int status;
	CPXENVptr env = CPXopenCPLEX(&status);
	if (status != 0) {
		throw std::runtime_error("Failed to open CPLEX environment: " + std::to_string(status));
	}
	for (const auto& inst:instances) {
		clock_t start = std::clock();
		std::cout << "Handling instance " << inst.first << std::endl;
		std::ifstream in("../instances/" + inst.first + ".tsp");
		if (!in) {
			std::cout << "Instance not found!" << std::endl;
			continue;
		}
		TSPInstance tsp(in);
		TSPSolution init = tspsolvers::solveGreedy(tsp);
		init = init.opt2();
		eliminateEdges(tsp, init.getCost(), env);
		std::cout << "Opt cost would be " << inst.second << std::endl;
		clock_t end = std::clock();
		double elapsed_secs = double(end-start)/CLOCKS_PER_SEC;
		std::cout << "Took " << elapsed_secs << " seconds" << std::endl;
	}
	CPXcloseCPLEX(&env);
}