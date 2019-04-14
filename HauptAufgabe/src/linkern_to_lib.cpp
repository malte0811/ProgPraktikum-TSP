#include <iostream>
#include <sstream>

int read() {
	std::string line;
	std::getline(std::cin, line);
	if (!std::cin) {
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

int main(int argc, char** argv) {
	if (argc!=2) {
		return 1;
	}
	std::cout << "NAME: " << argv[1] << std::endl
	<< "DIMENSION: " << read() << std::endl
	<< "TYPE: TOUR" << std::endl
	<< "TOUR_SECTION" << std::endl;
	int tmp = 0;
	while (tmp>=0) {
		tmp = read();
		if (tmp>=0) {
			std::cout << tmp+1 << std::endl;
		}
	}
	std::cout << "-1" << std::endl;
}