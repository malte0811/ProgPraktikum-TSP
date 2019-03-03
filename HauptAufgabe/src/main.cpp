#include <tsp_solvers.hpp>
#include <tsp_instance.hpp>
#include <fstream>
#include <map>
#include <string>

/**
 * Entfernt die Optionen (Parameter der Form --foo=bar) aus dem Vector der Argumente und gibt die Optionen als std::map
 * zurück
 * @param args Enthält vor dem Aufruf alle Argumente, nach dem Aufruf nur noch solche, die keine Optionen sind
 * @return die Optionen
 */
std::map<std::string, std::string> parseArgs(std::vector<std::string>& args) {
	std::map<std::string, std::string> ret;
	auto it = args.begin();
	while (it!=args.end()) {
		std::string arg = *it;
		bool erased = false;
		if (arg.size()>3 && arg[0]=='-' && arg[1]=='-') {
			arg = arg.substr(2);
			size_t equals = arg.find('=');
			if (equals!=std::string::npos && equals>0 && equals<arg.size()-1) {
				ret[arg.substr(0, equals)] = arg.substr(equals+1);
				it = args.erase(it);
				erased = true;
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
			throw std::runtime_error(retString+" is not a valid value for "+key);
		}
		return ret;
	} else {
		return defaultVal;
	}
}

int main(int argc, char** argv) {
	try {
		std::vector<std::string> args(argv+1, argv+argc);
		std::map<std::string, std::string> options = parseArgs(args);
		if (args.size()!=1 && args.size()!=2) {
			std::cerr << "Arguments: [options] <input file name> [<output file name>]" << std::endl;
			return 1;
		}
		bool useGreedy = getOption(options, "useGreedy", true);
		auto maxOpenSize = getOption<size_t>(options, "maxOpenSize", 1536*1024*1024);
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
		if (status!=0) {
			throw std::runtime_error("Failed to open CPLEX environment: "+std::to_string(status));
		}
		TSPInstance inst(in);
		in.close();
		TSPSolution initial;
		if (useGreedy) {
			initial = tspsolvers::solveGreedy(inst);
		}
		TSPSolution optimal = tspsolvers::solveLP(inst, useGreedy ? &initial : nullptr, env, maxOpenSize);
		if (args.size()==1) {
			optimal.write(std::cout);
		} else {
			std::ofstream out(argv[1]);
			if (!out) {
				std::cout << "Could not create/write to output file: " << argv[1] << std::endl;
				return 1;
			}
			optimal.write(out);
			out.close();
		}
		CPXcloseCPLEX(&env);
	} catch (std::runtime_error& err) {
		std::cout << "Error: " << err.what() << std::endl;
	}
}