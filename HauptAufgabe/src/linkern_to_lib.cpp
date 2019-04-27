#include <iostream>
#include <sstream>
#include <fstream>

int read(std::istream& in) {
	std::string line;
	std::getline(in, line);
	if (!in) {
		return -1;
	}
	std::stringstream ss(line);
	int ret;
	ss >> ret;
	if (!ss) {
		return -1;
	}
	return ret;
}

/*
 * Wandelt die Ausgabe von Concorde's linkern in eine TSPLIB-Tour um
 */
int main(int argc, char** argv) {
	if (argc != 3) {
		std::cout << "Needs exactly 2 arguments!" << std::endl;
		return 1;
	}
	std::ifstream in(argv[1]);
	std::ofstream out(argv[2]);
	if (!in || !out) {
		std::cout << "Input/Output files are invalid!" << std::endl;
		return 1;
	}
	out << "NAME: " << argv[2] << std::endl
		<< "DIMENSION: " << read(in) << std::endl
		<< "TYPE: TOUR" << std::endl
		<< "TOUR_SECTION" << std::endl;
	int tmp = 0;
	while (tmp >= 0) {
		tmp = read(in);
		if (tmp >= 0) {
			out << tmp + 1 << std::endl;
		}
	}
	out << "-1" << std::endl;
}