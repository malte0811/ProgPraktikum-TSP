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

using value_t = char;
using coeff_t = long;
using obj_t = typename std::common_type<value_t, coeff_t>::type;

class VariableRemover;

/**
 * Bestimmt die optimale Lösung eines ganzzahligen linearen Programms, das die folgenden Bedingungen erfüllt:
 * 1. Die Koeffizienten der Zielfunktion sind ganzzahlig
 * 2. Die Zielfunktion soll minimiert werden
 * 3. Alle Variablen sind binär, d.h. sollen in der Lösung den Wert 0 oder den Wert 1 annehmen
 */
class BranchAndCut {
public:
	BranchAndCut(LinearProgram& program, std::vector<CutGenerator *> gens, VariableRemover *remover, bool dfs);

	void setUpperBound(const std::vector<value_t>& value, obj_t cost);

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
		explicit SystemBounds(variable_id varCount);

		SystemBounds(const BranchAndCut::SystemBounds& old, const std::vector<variable_id>& idMap,
					 variable_id newVarCount);

		VariableBounds operator[](variable_id id) const;

		void fix(variable_id var, LinearProgram::BoundType b);

		void setBound(variable_id var, LinearProgram::BoundType b, value_t newBound);

		bool operator==(const SystemBounds& other) const;

		variable_id getVarCount() const;

	private:
		std::vector<bool> fixUpper;
		std::vector<bool> fixLower;
	};

	struct BranchNode {
		BranchNode() = delete;

		SystemBounds bounds;
		double value;

		bool operator<(const BranchNode& other) const;

		bool operator==(const BranchNode& other) const;
	};

	LinearProgram& problem;
	//Die aktuelle obere Schranke
	obj_t upperBound;
	//Die Anzahl der Variablen im LP
	variable_id varCount;
	//Die Variablenbelegung, mit der die obere Schranke erreicht wird
	std::vector<value_t> currBest;
	//Die aktuelle fraktionale Lösung
	LinearProgram::Solution fractOpt;
	//Die Einträge geben an, wie lange eine gegebenen Ungleichung nicht mehr mit Gleichheit erfüllt war
	std::vector<size_t> sinceSlack0;
	//Die Koeffizienten der Zielfunktion
	std::vector<coeff_t> objCoefficients;
	//Die Menge der offenen Knoten, sortiert nach dem Wert der Zielfunktion
	std::multiset<BranchNode> open;
	//Die aktuellen Variablenschranken
	std::vector<VariableBounds> currentBounds;

	std::vector<LinearProgram::Constraint> recentlyRemoved;
	//Gibt an, ob der Suchbaum nach der üblichen Strategie oder mit DFS durchlaufen werden soll
	const bool dfs;
	//Die Schnittebenen-Generatoren, die erfüllt sein müssen
	const std::vector<CutGenerator *> generators;
	VariableRemover *remover;
	//Anzahl der Nebenbedingungen am Start der Berechnung
	const size_t constraintsAtStart;
	//Toleranz für allgemeine Verwendung
	const RelativeTolerance generalTolerance;
	//Toleranz, nach der entschieden wird, ob ein Wert ganzzahlig ist
	const lemon::Tolerance<double> intTolerance;
	//Anzahl der bereits behandelten Knoten im "Suchbaum"
	size_t handledNodes = 0;
	//Die Knoten im Suchbaum, die aktuell behandelt werden
	std::vector<BranchNode *> currStack;
	/*
	 * Ordnet jeder Variablen ihre ursprüngliche ID zu. Wird genutzt, um Variablen über Entfernungsschritte hinweg
	 * eindeutig identifizieren zu können
	 */
	std::vector<variable_id> correspondingOrigVar;

	bool canImproveBound(double fractBound);

	void branchAndBound(BranchNode& node, bool isRoot);

	void branch(variable_id variable, value_t val, LinearProgram::BoundType bound, const SystemBounds& parent,
				double objValue, bool immediate);

	void solveLP(LinearProgram::Solution& out, bool isRoot);

	void countSolutionSlack(const LinearProgram::Solution& sol);

	void setupBounds(const SystemBounds& bounds);

	CutGenerator::CutStatus readdRemovedConstraints(const LinearProgram::Solution& sol);

	bool cleanupOldConstraints();

	void removeVariables(const std::vector<variable_id>& toRemove);
};


#endif