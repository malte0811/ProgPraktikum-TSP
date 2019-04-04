#include <tsp_instance.hpp>
#include <sstream>
#include <cmath>
#include <array>
#include <cassert>

/**
 * Liest einen Wert aus dem gegebenen istrean und wirft einen Fehler, falls dies nicht möglich ist
 * @tparam T Der Typ des zu lesenden Wertes
 */
template<typename T>
T readOrThrow(std::istream& input) {
	T ret;
	input >> ret;
	if (!input) {
		throw std::invalid_argument("Could not read input!");
	}
	return ret;
}

/**
 * Liest eine STSP im TSPLIB-Format ein
 * @param input Die Quelle der Eingabe
 */
TSPInstance::TSPInstance(std::istream& input) {
	std::string line;
	EdgeWeightType edgeType = euc_2d;
	EdgeFormat edgeFormat = full_matrix;
	bool emptyLines = false;
	while (std::getline(input, line)) {
		if (!line.empty()) {
			if (emptyLines) {
				std::cout << "Skipped empty line(s)" << std::endl;
			}
			std::stringstream ss(line);
			std::string keyword;
			/*
			 * Doppelpunkt ignorieren. Der Doppelpunkt kann sowohl direkt hinter dem "Operator" als auch durch Whitespace
			 * abgetrennt auftreten.
			 */
			ss >> keyword >> std::ws;
			if (ss.peek()==':') {
				ss.ignore();
			} else if (keyword.back()==':') {
				keyword = keyword.substr(0, keyword.size()-1);
			}
			if (keyword=="NAME") {
				ss >> name;
			} else if (keyword=="COMMENT") {
				//NOP
			} else if (keyword=="DIMENSION") {
				auto nodeCount = readOrThrow<city_id>(ss);
				distances.resize(static_cast<size_t>(nodeCount-1));
				for (size_t i = 0; i<static_cast<size_t>(nodeCount-1); ++i) {
					distances[i].resize(i+1);
				}
			} else if (keyword=="EDGE_WEIGHT_TYPE") {
				auto type = readOrThrow<std::string>(ss);
				if (type=="EUC_2D") {
					edgeType = euc_2d;
				} else if (type=="CEIL_2D") {
					edgeType = ceil_2d;
				} else if (type=="EXPLICIT") {
					edgeType = explicit_;
				} else if (type=="GEO") {
					edgeType = geo;
				} else if (type=="ATT") {
					edgeType = att;
				} else {
					throw std::runtime_error("Unknown edge weight type: "+type);
				}
			} else if (keyword=="EDGE_WEIGHT_FORMAT") {
				auto type = readOrThrow<std::string>(ss);
				if (type=="FULL_MATRIX") {
					edgeFormat = full_matrix;
				} else if (type=="LOWER_DIAG_ROW") {
					edgeFormat = lower_diag_row;
				} else if (type=="UPPER_DIAG_ROW") {
					edgeFormat = upper_diag_row;
				} else if (type=="UPPER_ROW") {
					edgeFormat = upper_row;
				} else if (type!=
						   "FUNCTION") {//Function kommt in einigen Instanzen vor, die sonst keine ungewöhnlichen Daten haben (gr431)
					throw std::runtime_error("Unknown edge format: "+type);
				}
			} else if (keyword=="NODE_COORD_SECTION") {
				readNodes(input, edgeType);
				input >> std::ws;
			} else if (keyword=="EDGE_WEIGHT_SECTION") {
				readEdges(input, edgeFormat);
				input >> std::ws;
			} else if (keyword=="TYPE") {
				auto type = readOrThrow<std::string>(ss);
				if (type!="TSP") {
					throw std::invalid_argument("Input is not a symmetric TSP instance!");
				}
			} else if (keyword=="EOF") {
				break;
			} else {
				std::cout << "Unknown keyword found, ignoring: " << keyword << std::endl;
			}
		} else {
			emptyLines = true;
		}
	}
}

/**
 * Wandelt einen Wert der Form ggg.mm (ggg Grd und mm Minuten) in Bogenmaß um
 */
double coordToLatLong(double val) {
	//TSPLIB-Beschreibung nutzt hier nint (round), das geht aber nicht für >=50 Minuten
	int deg = static_cast<int>(val);
	double min = val-deg;
	return M_PI*(deg+5.0*min/3.0)/180.0;
}

void TSPInstance::readNodes(std::istream& input, EdgeWeightType type) {
	//Knotenkoordinaten einlesen (Es werden nur Formate mit 2 Koordinaten pro Knoten unterstützt)
	struct Vec2 {
		double x, y;
	};
	city_id nodeCount = getCityCount();
	std::vector<Vec2> nodeLocations(static_cast<size_t>(nodeCount));
	//Wurde die Position eines gegebenen Knotens gesetzt?
	std::vector<bool> set(static_cast<size_t>(nodeCount), false);
	city_id setCount = 0;
	for (city_id iteration = 0; iteration<nodeCount; ++iteration) {
		auto id = readOrThrow<city_id>(input)-1;
		assert(id<nodeCount);
		nodeLocations[id] = {readOrThrow<double>(input), readOrThrow<double>(input)};
		if (set[id]) {
			throw std::runtime_error("Location for city "+std::to_string(id)+" was set twice");
		}
		set[id] = true;
		++setCount;
	}
	if (setCount<nodeCount) {
		throw std::runtime_error(
				"Only "+std::to_string(setCount)+" of "+std::to_string(nodeCount)+" node locations have been set");
	}
	//Über alle Kanten iterieren und die Distanzen berechnen
	for (city_id higherId = 1; higherId<nodeCount; ++higherId) {
		for (city_id lowerId = 0; lowerId<higherId; ++lowerId) {
			cost_t distance = 0;
			double distX = nodeLocations[lowerId].x-nodeLocations[higherId].x;
			double distY = nodeLocations[lowerId].y-nodeLocations[higherId].y;
			switch (type) {
				case euc_2d:
					distance = static_cast<cost_t>(std::round(std::sqrt(distX*distX+distY*distY)));
					break;
				case ceil_2d:
					distance = static_cast<cost_t>(std::ceil(std::sqrt(distX*distX+distY*distY)));
					break;
				case explicit_:
					throw std::runtime_error("Weight type is EXPLICIT, but a NODE_COORD_SECTION exists!");
				case geo: {
					double lat1 = coordToLatLong(nodeLocations[lowerId].x);
					double long1 = coordToLatLong(nodeLocations[lowerId].y);
					double lat2 = coordToLatLong(nodeLocations[higherId].x);
					double long2 = coordToLatLong(nodeLocations[higherId].y);
					//Aus der TSPLib-Dokumentation
					const double RRR = 6378.388;
					double q1 = std::cos(long1-long2);
					double q2 = std::cos(lat1-lat2);
					double q3 = std::cos(lat1+lat2);
					distance = static_cast<cost_t>(RRR*std::acos(0.5*((1.0+q1)*q2-(1.0-q1)*q3))+1.0);
				}
					break;
				case att: {
					//Aus der TSPLib-Dokumentation
					double rij = sqrt((distX*distX+distY*distY)/10.0);
					auto tij = static_cast<cost_t>(std::lround(rij));
					if (tij<rij) {
						distance = tij+1;
					} else {
						distance = tij;
					}
				}
					break;
			}
			distances[higherId-1][lowerId] = distance;
		}
	}
}

void TSPInstance::readEdges(std::istream& input, TSPInstance::EdgeFormat type) {
	const city_id nodeCount = getCityCount();
	city_id rowCount = nodeCount;
	if (type==upper_row) {
		rowCount--;
	}
	for (city_id row = 0; row<rowCount; ++row) {
		//Die Spalte, der die erste Zahl entspricht
		city_id minCol = 0;
		//Die Anzahl der Spalten in dieser Zeile
		city_id colCount = 0;
		switch (type) {
			case full_matrix:
				colCount = nodeCount;
				break;
			case lower_diag_row:
				colCount = row+1;
				break;
			case upper_diag_row:
				minCol = row;
				colCount = nodeCount-row;
				break;
			case upper_row:
				minCol = row+1;
				colCount = nodeCount-row-1;
				break;
		}
		for (city_id col = minCol; col<minCol+colCount; ++col) {
			setDistance(row, col, readOrThrow<cost_t>(input));
		}
	}
}

void TSPInstance::setDistance(city_id a, city_id b, cost_t dist) {
	if (a>b) {
		distances[a-1][b] = dist;
	} else if (b>a) {
		distances[b-1][a] = dist;
	}
}

std::string TSPInstance::getName() const {
	return name;
}
