#ifndef SUBTOUR_CUT_GEN_HPP
#define SUBTOUR_CUT_GEN_HPP

#include <cut_generator.hpp>
#include <lemon/smart_graph.h>
#include <lemon_fixes/nagamochi_ibaraki.h>
#include <tsp_instance.hpp>

class SubtourCutGen : public CutGenerator {
public:
	explicit SubtourCutGen(const TSPInstance& inst);

	CutStatus validate(LinearProgram& lp, const std::vector<double>& solution) override;

private:
	void addConnectivityConstraints(LinearProgram& lp);

	void addCutConstraint(LinearProgram& lp);

	//Die TSP-Instanz
	const TSPInstance& tsp;
	//Der Arbeitsgraph. Enthält alle Kanten, deren Wert in der LP-Lösung nicht 0 ist
	Graph workGraph;
	//Der Grundzustand des Arbeitsgraphen: Keine Kanten, nur die korrekte Anzahl an Knoten
	Graph::Snapshot baseState;
	//Ordnet jedem Knoten im TSP-Graphen einen Knoten im Arbeitsgraphen zu
	Graph::NodeMap <Graph::Node> origToWork;
	//Die Umkehrabbildung zu origToWork
	Graph::NodeMap <Graph::Node> workToOrig;
	//Die Kantenkapazitäten im Arbeitsgraphen
	Graph::EdgeMap<double> capacity;
	//Der Min-Cut-Löser
	lemon::NagamochiIbaraki<Graph, Graph::EdgeMap<double>> minCut;
	lemon::Tolerance<double> tolerance;
};


#endif