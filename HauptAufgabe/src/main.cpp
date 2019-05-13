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

using std::size_t;

std::vector<std::string> splitOnChar(const std::string& string, char delim) {
	std::vector<std::string> ret;
	std::stringstream in(string);
	std::string segment;
	while (std::getline(in, segment, delim)) {
		ret.push_back(segment);
	}
	return ret;
}//TODO: Warum verwendest Du diese Funktion nicht in parseArgs?

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
			size_t equals = arg.find('=');
			if (equals != std::string::npos && equals > 0 && equals < arg.size() - 1) {
				ret[arg.substr(0, equals)] = arg.substr(equals + 1);
				it = args.erase(it);
				erased = true;
			}
		}
		if (!erased) {
			++it;
		}
	}
	return ret;
}//TODO: Wenn es wirklich fancy sein soll dann nimmst Du iteratoren als argument
// und shiftest die zu loeschenden elemente ans ende, Zum schluss gibst Du dann einen Iterator auf das erste ungueltige element
//zurueck. Dann ist parseArgs unabhaengig von der Datenstruktur, in der die strings gespeichert sind.

//TODO: Mir ist unklar, warum hier die Optionen geloescht werden.
// Die koennen ja auch einfach am Ende mit der Map geloescht werden.
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
//TODO: Genereller Tipp: Leerzeilen zwischen zusammenhaengenden Codebloecken verbessern die Lesbarkeit.
int main(int argc, char **argv) {
	try {
		std::vector<std::string> args(argv + 1, argv + argc);
		std::map<std::string, std::string> options = parseArgs(args);
		if (args.size() != 1 && args.size() != 2) {
		   //TODO: Btw. in C++ kann man and or und not verwenden (alias fuer &&, ||, !). Das ist einfacher zu lesen.
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
		int status;
		CPXENVptr env = CPXopenCPLEX(&status);
		if (status != 0) {
			throw std::runtime_error("Failed to open CPLEX environment: " + std::to_string(status));
		}//TODO: Warum ist das hier so eingeschoben. Lieber vorm oeffnen des Files oder nach dem Einlesen der Instanz.
		TSPInstance inst(in);
		in.close();
		TSPSolution initial;
		if (startingBound == "<greedy>" || startingBound == "<greedy2opt>") {
			initial = tspsolvers::solveGreedy(inst);
			if (startingBound == "<greedy2opt>") {
				initial = initial.opt2();
			}
		} else if (startingBound != noBound) {
			std::ifstream boundIn(startingBound);
			initial = TSPSolution(inst, boundIn);
		}
		TSPSolution optimal = tspsolvers::solveLP(inst, initial.isValid() ? &initial : nullptr, env, useDFS,
												  generators);
		if (args.size() == 1) {
			optimal.write(std::cout);
		} else {
			std::ofstream out(args[1]);
			if (!out) {
				std::cout << "Could not create/write to output file: " << args[1] << std::endl;
            //TODO: Hier ist es evtl sinnvoll das Ergebnis nach std::cout zu schreiben.
            // Sonst rechnet man evtl ewig und dann bekommt man keinen Output, weil der den File nicht oeffnen konnte.
            return 1;
         }
         optimal.write(out);
         out.close();
		}
		if (expectedValue > 0 && optimal.getCost() != expectedValue) {
			std::cerr << "Found tour of cost " << optimal.getCost() << ", but expected cost " << expectedValue
					  << std::endl;
		}
		CPXcloseCPLEX(&env);
	} catch (std::runtime_error& err) {
		std::cout << "Error: " << err.what() << std::endl;
	}
}
