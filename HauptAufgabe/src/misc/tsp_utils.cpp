#include <tsp_utils.hpp>
#include <lemon/core.h>
#include <cstddef>
#include <iostream>
#include <tsp_lp_data.hpp>
#include <linear_program.hpp>
#include <tsp_instance.hpp>

/**
 * Fügt die Kanten, deren Variablen echt positiv sind, zum Graphen hinzu.
 * @param lpData Gibt an, welche Kanten welchen Variablen entsprechen
 * @param tolerance Gibt an, wann eine Variable als positiv gilt
 * @param solution Die LP-Lösung, zu der der Support-Graph erstellt werden soll
 * @param workGraph Der Ausgabegraph
 * @param origToWork Ordnet der Städten in der TSP-Instanz Knoten zu
 * @param workToOrig Ordnet den Knoten Städte in der TSP-Instanz zu
 * @param c Die Kosten der Kanten bzw. die Werte der entsprechenden Variablen
 */
void tsp_util::addSupportGraphEdges(const TspLpData& lpData, lemon::Tolerance<double> tolerance,
									const std::vector<double>& solution, Graph& workGraph,
									std::vector<Graph::Node>& origToWork, Graph::NodeMap <city_id>& workToOrig,
									Graph::EdgeMap<double>& c) {
	if (Graph::NodeIt(workGraph) == lemon::INVALID) {
		//Falls der Graph noch keine Knoten enthält, wird die richtige Menge an Knoten hinzugefügt
		origToWork.resize(lpData.getTSP().getCityCount());
		for (city_id i = 0; i < lpData.getTSP().getCityCount(); ++i) {
			Graph::Node newNode = workGraph.addNode();
			workToOrig[newNode] = i;
			origToWork[i] = newNode;
		}
	} else {
		//Bereits existierende Kanten entfernen
		for (Graph::EdgeIt it(workGraph); it != lemon::INVALID;) {
			Graph::Edge e = it;
			++it;
			workGraph.erase(e);
		}
	}
	/*
	 * Alle Kanten einfügen, deren Variablen einen echt positiven Wert haben
	 */
	for (variable_id i = 0; i < lpData.getVariableCount(); ++i) {
		if (tolerance.positive(solution[i])) {
			const TspLpData::Edge& e = lpData.getEdge(i);
			Graph::Edge inWork = workGraph.addEdge(origToWork[e.first], origToWork[e.second]);
			c[inWork] = solution[i];
		}
	}
}

std::string tsp_util::readKeyword(std::istream& in) {
	std::string keyword;
	in >> keyword >> std::ws;
	/*
	 * Doppelpunkt ignorieren. Der Doppelpunkt kann sowohl direkt hinter dem "Operator" als auch durch Whitespace
	 * abgetrennt auftreten.
	 */
	if (in.peek() == ':') {
		in.ignore();
	} else if (keyword.back() == ':') {
		keyword = keyword.substr(0, keyword.size() - 1);
	}
	return keyword;
}

/**
 * Vergleicht Constraints, bei denen die Nonzeroes in aufsteigender ID-Reihenfolge angegeben sind. Constraints mit
 * kleinerem Slack sind dabei immer klener.
 */
bool tsp_util::CompareOrderedConstraint::operator()(const ConstraintWithSlack& p1, const ConstraintWithSlack& p2) {
	if (p1.slack != p2.slack) {
		return p1.slack < p2.slack;
	}
	const LinearProgram::Constraint& c1 = p1.constraint;
	const LinearProgram::Constraint& c2 = p2.constraint;
	if (c1.getNonzeroes().size() != c2.getNonzeroes().size()) {
		return c1.getNonzeroes().size() < c2.getNonzeroes().size();
	}
	if (c1.getRHS() != c2.getRHS()) {
		return c1.getRHS() < c2.getRHS();
	}
	if (c1.getSense() != c2.getSense()) {
		return c1.getSense() < c2.getSense();
	}
	for (size_t i = 0; i < c1.getNonzeroes().size(); ++i) {
		if (c1.getNonzeroes()[i] != c2.getNonzeroes()[i]) {
			return c1.getNonzeroes()[i] < c2.getNonzeroes()[i];
		}
		if (c1.getCoeffs()[i] != c2.getCoeffs()[i]) {
			return c1.getCoeffs()[i] < c2.getCoeffs()[i];
		}
	}
	return false;
}
