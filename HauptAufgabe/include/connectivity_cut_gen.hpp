#ifndef CONNECTIVITY_CUT_GEN_HPP
#define CONNECTIVITY_CUT_GEN_HPP

#include <lemon/tolerance.h>
#include <cut_generator.hpp>
#include <tsp_instance.hpp>
#include <vector>

class LinearProgram;

class TspLpData;

/**
 * Fügt in nicht zusammenhängenden Graphen die Subtour-Constraints zu den Zusammenhangskomponenten hinzu
 */
class ConnectivityCutGen : public CutGenerator {
public:
	explicit ConnectivityCutGen(const TSPInstance& inst, const TspLpData& lpData);

	CutStatus validate(LinearProgram& lp, const std::vector<double>& solution, CutStatus) override;

private:
	//Der Arbeitsgraph. Enthält alle Kanten, deren Wert in der LP-Lösung nicht 0 ist
	Graph workGraph;
	//Ordnet jedem Stadt der TSP-Instanz einen Knoten im Arbeitsgraphen zu
	std::vector<Graph::Node> origToWork;
	//Die Umkehrabbildung zu origToWork
	Graph::NodeMap <city_id> workToOrig;
	const TSPInstance& tsp;
	const TspLpData& lpData;
	lemon::Tolerance<double> tolerance;
};

#endif