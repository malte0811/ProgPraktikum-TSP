
#include <tsp_instance.hpp>
#include <comb_cut_gen.hpp>

int main() {
	std::stringstream in("NAME: Test\n"
						 "TYPE: TSP\n"
						 "DIMENSION: 12\n"
						 "EDGE_WEIGHT_FORMAT: LOWER_DIAG_ROW\n"
						 "EDGE_WEIGHT_SECTION\n"
						 "1\n"
						 "1 1\n"
						 "1 1 1\n"
						 "1 1 1 1\n"
						 "1 1 1 1 1\n"
						 "1 1 1 1 1 1\n"
						 "1 1 1 1 1 1 1\n"
						 "1 1 1 1 1 1 1 1\n"
						 "1 1 1 1 1 1 1 1 1\n"
						 "1 1 1 1 1 1 1 1 1 1\n"
						 "1 1 1 1 1 1 1 1 1 1 1\n"
						 "1 1 1 1 1 1 1 1 1 1 1 1\n"
	);
	TSPInstance tsp(in);
	std::vector<double> sol(tsp.getEdgeCount(), 0);
	sol[tsp.getVariable(0, 1)] = 0.75;
	sol[tsp.getVariable(4, 0)] = 0.75;
	sol[tsp.getVariable(5, 0)] = 0.5;
	sol[tsp.getVariable(4, 5)] = 0.25;
	sol[tsp.getVariable(1, 5)] = 0.25;

	sol[tsp.getVariable(6, 5)] = 0.75;
	sol[tsp.getVariable(6, 7)] = 0.25;
	sol[tsp.getVariable(2, 5)] = 0.25;
	sol[tsp.getVariable(2, 7)] = 0.75;

	sol[tsp.getVariable(8, 7)] = 0.5;
	sol[tsp.getVariable(3, 7)] = 0.5;
	sol[tsp.getVariable(3, 8)] = 0.5;

	//Odd nodes
	sol[tsp.getVariable(1, 9)] = 1;
	sol[tsp.getVariable(2, 9)] = 1;
	sol[tsp.getVariable(3, 10)] = 1;
	sol[tsp.getVariable(4, 10)] = 1;
	sol[tsp.getVariable(6, 11)] = 1;
	sol[tsp.getVariable(8, 11)] = 1;
	int status;
	CPXENVptr env = CPXopenCPLEX(&status);
	LinearProgram lp(env, "temp", LinearProgram::minimize);
	tsp.setupBasicLP(lp);
	CombCutGen combGen(tsp);
	std::cout << combGen.validate(lp, sol, CutGenerator::valid) << std::endl;
	CPXcloseCPLEX(&env);
}