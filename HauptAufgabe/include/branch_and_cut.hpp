#ifndef BRANCH_AND_CUT_HPP
#define BRANCH_AND_CUT_HPP


#include <lemon/tolerance.h>
#include <linear_program.hpp>
#include <cut_generator.hpp>
#include <queue>
#include <stack>
#include <iostream>
#include <array>
#include <map>
#include <set>

/**
 * Bestimmt die optimale Lösung eines ganzzahligen linearen Programms, bei dem die Zielfunktion ganzzahlige
 * Koeffizienten hat
 */
class BranchAndCut {
public:
	BranchAndCut(LinearProgram& program, const std::vector<CutGenerator*>& gens, size_t maxOpenSize);

	void setUpperBound(const std::vector<long>& value, long cost);

	std::vector<long> solve();

private:

	/**
	 * Speichert die Beschränkungen für eine einzelne LP-Variable
	 */
	struct VariableBounds {
		long min, max;

		long& operator[](LinearProgram::BoundType b);
	};

	/**
	 * Speichert die Beschränkungen für alle Variablen des LPs
	 */
	struct SystemBounds {
		std::map<variable_id, VariableBounds> bounds;
		std::vector<bool> fixUpper;
		std::vector<bool> fixLower;
		BranchAndCut* owner;

		explicit SystemBounds(BranchAndCut* owner);

		VariableBounds operator[](variable_id id) const;

		void fix(variable_id var, LinearProgram::BoundType b);
	};

	struct BranchNode {
		BranchNode() = delete;

		SystemBounds bounds;
		double value;
		LinearProgram::Goal goal;

		bool operator<(const BranchNode& other) const;

		size_t estimateSize() const;
	};

	LinearProgram& problem;
	//Die aktuelle obere Schranke
	double upperBound;
	//Die Anzahl der Variablen im LP
	const variable_id varCount;
	//Das Ziel des LP's
	const LinearProgram::Goal goal;
	//Die Variablenbelegung, mit der die obere Schranke erreicht wird
	std::vector<long> currBest;
	//Die aktuelle fraktionale Lösung
	LinearProgram::Solution fractOpt;
	//Die Einträge geben an, wie lange eine gegebenen Ungleichung nicht mehr mit Gleichheit erfüllt war
	std::vector<size_t> sinceSlack0;
	//Die Koeffizienten der Zielfunktion
	std::vector<long> objCoefficients;
	//Die Menge der offenen Knoten, sortiert nach dem Wert der Zielfunktion
	std::set<BranchNode> open;
	//Die Variablenschranken im ursprünglichen LP
	std::vector<VariableBounds> defaultBounds;
	//Die aktuellen Variablenschranken
	std::vector<VariableBounds> currentBounds;
	//Die geschätzte Größe der offenen Menge in Bytes
	size_t openSize = 0;
	//Die maximale Größe der offenen Menge
	const size_t maxOpenSize;
	//Die Schnittebenen-Generatoren, die erfüllt sein müssen
	const std::vector<CutGenerator*> generators;
	//Anzahl der Nebenbedingungen am Start der Berechnung
	const size_t constraintsAtStart;
	//Toleranz für allgemeine Verwendung
	const lemon::Tolerance<double> generalTolerance;
	//Toleranz, nach der entschieden wird, ob ein Wert ganzzahlig ist
	const lemon::Tolerance<double> intTolerance;

	static bool isBetter(double a, double b, LinearProgram::Goal goal);

	void branchAndBound(BranchNode& node, bool dfs);

	void branch(variable_id variable, long val, LinearProgram::BoundType bound,
				const SystemBounds& parent, double objValue, bool immediate, bool dfs);

	void solveLP(LinearProgram::Solution& out);

	void countSolutionSlack(const LinearProgram::Solution& sol);

	void setupBounds(const SystemBounds& bounds);
};


#endif