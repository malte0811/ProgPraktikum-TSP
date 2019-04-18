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
#include <relative_tolerance.hpp>

using value_t = long;
using coeff_t = long;

class VariableRemover;
/**
 * Bestimmt die optimale Lösung eines ganzzahligen linearen Programms, bei dem die Zielfunktion ganzzahlige
 * Koeffizienten hat
 * TODO evtl auf 0-1-Programme einschränken? Dann ist die map im VariableBounds nicht mehr nötig
 */
class BranchAndCut {
public:
	BranchAndCut(LinearProgram& program, std::vector<CutGenerator *>  gens, VariableRemover *remover,
				 size_t maxOpenSize);

	void setUpperBound(const std::vector<value_t>& value, value_t cost);

	std::vector<value_t> solve();

private:

	/**
	 * Speichert die Beschränkungen für eine einzelne LP-Variable
	 */
	struct VariableBounds {
		value_t min, max;

		value_t& operator[](LinearProgram::BoundType b);

		bool operator==(const VariableBounds& other) const;

		bool isFixed() const;
	};

	/**
	 * Speichert die Beschränkungen für alle Variablen des LPs
	 */
	class SystemBounds {
	public:
		explicit SystemBounds(BranchAndCut* owner);

		SystemBounds(const BranchAndCut::SystemBounds& old,
					 const std::vector<variable_id>& idMap, variable_id newVarCount);

		VariableBounds operator[](variable_id id) const;

		void fix(variable_id var, LinearProgram::BoundType b);

		void setBound(variable_id var, LinearProgram::BoundType b, value_t newBound);

		bool operator==(const SystemBounds& other) const;

		size_t getFixedCount() const;

		size_t estimateSize() const;

	private:
		std::map<variable_id, VariableBounds> bounds;
		std::vector<bool> fixUpper;
		std::vector<bool> fixLower;
		BranchAndCut* owner;
		size_t fixedCount;
	};

	struct BranchNode {
		BranchNode() = delete;

		SystemBounds bounds;
		double value;
		LinearProgram::Goal goal;

		bool operator<(const BranchNode& other) const;

		size_t estimateSize() const;

		size_t getFixedCount() const;

		bool operator==(const BranchNode& other) const;
	};

	LinearProgram& problem;
	//Die aktuelle obere Schranke
	double upperBound;
	//Die Anzahl der Variablen im LP
	variable_id varCount;
	//Das Ziel des LP's
	const LinearProgram::Goal goal;
	//Die Variablenbelegung, mit der die obere Schranke erreicht wird
	std::vector<value_t> currBest;
	//Die aktuelle fraktionale Lösung
	LinearProgram::Solution fractOpt;
	//Die Einträge geben an, wie lange eine gegebenen Ungleichung nicht mehr mit Gleichheit erfüllt war
	std::vector<size_t> sinceSlack0;
	//Die Koeffizienten der Zielfunktion
	std::vector<coeff_t> objCoefficients;
	//Die Menge der offenen Knoten, sortiert nach dem Wert der Zielfunktion
	//TODO für max. muss anders sortiert werden. Will ich max überhaupt zulassen?
	std::multiset<BranchNode> open;
	//Die Variablenschranken im ursprünglichen LP
	std::vector<VariableBounds> defaultBounds;
	//Die aktuellen Variablenschranken
	std::vector<VariableBounds> currentBounds;
	//Die geschätzte Größe der offenen Menge in Bytes
	size_t openSize = 0;

	std::vector<LinearProgram::Constraint> recentlyRemoved;
	//Die maximale Größe der offenen Menge
	const size_t maxOpenSize;
	//Die Schnittebenen-Generatoren, die erfüllt sein müssen
	const std::vector<CutGenerator*> generators;
	VariableRemover *remover;
	//Anzahl der Nebenbedingungen am Start der Berechnung
	const size_t constraintsAtStart;
	//Toleranz für allgemeine Verwendung
	const RelativeTolerance generalTolerance;
	//Toleranz, nach der entschieden wird, ob ein Wert ganzzahlig ist
	const lemon::Tolerance<double> intTolerance;
	//Anzahl der bereits behandelten Knoten im "Suchbaum"
	size_t handledNodes = 0;
	BranchNode* currentNode = nullptr;

	bool isBetter(double a, double b);

	void branchAndBound(BranchNode& node, bool dfs, bool isRoot);

	void branch(variable_id variable, value_t val, LinearProgram::BoundType bound,
				const SystemBounds& parent, double objValue, bool immediate, bool dfs);

	void solveLP(LinearProgram::Solution& out, bool isRoot);

	void countSolutionSlack(const LinearProgram::Solution& sol);

	void setupBounds(const SystemBounds& bounds);

	CutGenerator::CutStatus readdRemovedConstraints(const LinearProgram::Solution& sol);

	bool cleanupOldConstraints();

	void removeVariables(const std::vector<variable_id>& toRemove);
};


#endif