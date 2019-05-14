#include <cstddef>
#include <iostream>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <tsp_instance.hpp>
#include <tsp_solvers.hpp>
#include <utility>
#include <vector>
#include <ilcplex/cplex.h>
#include <ilcplex/cpxconst.h>
#include <tsp_solution.hpp>
#include <linear_program.hpp>

using std::size_t;

std::vector<std::string> splitOnChar(const std::string& string, char delim) {
	std::vector<std::string> ret;
	std::stringstream in(string);
	std::string segment;
	while (std::getline(in, segment, delim)) {
		ret.push_back(segment);
	}
	return ret;
}

/**
 * Entfernt die Optionen (Parameter der Form --foo=bar) aus dem Vector der Argumente und gibt die Optionen als std::map
 * zurück
 * @param args Enthält vor dem Aufruf alle Argumente, nach dem Aufruf nur noch solche, die keine Optionen sind
 * @return die Optionen
 */
std::map<std::string, std::string> parseArgs(std::vector<std::string>& args) {
	std::map<std::string, std::string> ret;
	auto it = args.begin();
	while (it != args.end()) {
		std::string arg = *it;
		bool erased = false;
		if (arg.size() > 3 && arg[0] == '-' && arg[1] == '-') {
			arg = arg.substr(2);
			std::vector<std::string> split = splitOnChar(arg, '=');
			if (split.size() == 2) {
				ret[split[0]] = split[1];
				it = args.erase(it);
				erased = true;
			} else {
				throw std::runtime_error("Invalid argument: --" + arg);
			}
		}
		if (!erased) {
			++it;
		}
	}
	return ret;
}

/**
 * Liest den Wert der angegebenen Option aus, entfernt ihn aus der map und gibt ihn zurück
 * @tparam T Der Typ des Wertes der Option
 * @param options Alle Optionen
 * @param key Der Name der Option
 * @param defaultVal Der Standardwert (falls die Option nicht angegeben wurde)
 * @return Den Wert der Option
 */
template<typename T>
T getOption(std::map<std::string, std::string>& options, const std::string& key, T defaultVal) {
	if (options.count(key)) {
		std::string retString = options[key];
		std::stringstream valStream(retString);
		options.erase(key);
		T ret;
		valStream >> ret;
		if (!valStream) {
			throw std::runtime_error(retString + " is not a valid value for " + key);
		}
		return ret;
	} else {
		return defaultVal;
	}
}

int main(int argc, char **argv) {
	try {
		std::vector<std::string> args(argv + 1, argv + argc);
		std::map<std::string, std::string> options = parseArgs(args);
		if (args.size() != 1 && args.size() != 2) {
			std::cerr << "Arguments: [options] <input file name> [<output file name>]" << std::endl;
			return 1;
		}

		const std::string noBound = "<none>";
		std::string startingBound = getOption<std::string>(options, "startingBound", "<greedy>");
		std::string generatorString = getOption<std::string>(options, "cutGens", tspsolvers::cutgens::defaultGens);
		std::vector<std::string> generators = splitOnChar(generatorString, ',');
		bool useDFS = getOption(options, "dfs", false);
		cost_t expectedValue = getOption(options, "expectedResult", 0U);
		if (!options.empty()) {
			std::cerr << "Found unknown options:" << std::endl;
			for (const auto& entry:options) {
				std::cerr << "--" << entry.first << "=" << entry.second << std::endl;
			}
			return 1;
		}

		std::ifstream in(args[0]);
		if (!in) {
			std::cout << "File does not exist: " << args[0] << std::endl;
			return 1;
		}
		TSPInstance inst(in);
		in.close();

		TSPSolution initial;
		if (startingBound == "<greedy>" || startingBound == "<greedy2opt>") {
			initial = tspsolvers::solveGreedy(inst);
			if (startingBound == "<greedy2opt>") {
				initial = tspsolvers::opt2(initial, inst);
			}
		} else if (startingBound != noBound) {
			std::ifstream boundIn(startingBound);
			initial = TSPSolution(inst, boundIn);
		}

		SharedCplexEnv env = LinearProgram::openCPLEX();
		TSPSolution optimal = tspsolvers::solveLP(inst, initial.isValid() ? &initial : nullptr, env, useDFS,
												  generators);

		bool outputToConsole = args.size() == 1;
		if (!outputToConsole) {
			std::ofstream out(args[1]);
			if (out) {
				optimal.write(out);
				out.close();
			} else {
				std::cout << "Could not create/write to output file: " << args[1] << ", printing to console instead"
						  << std::endl;
				outputToConsole = true;
			}
		}
		if (outputToConsole) {
			optimal.write(std::cout);
		}
		if (expectedValue > 0 && optimal.getCost() != expectedValue) {
			std::cerr << "Found tour of cost " << optimal.getCost() << ", but expected cost " << expectedValue
					  << std::endl;
		}
	} catch (std::runtime_error& err) {
		std::cout << "Error: " << err.what() << std::endl;
	}
}
