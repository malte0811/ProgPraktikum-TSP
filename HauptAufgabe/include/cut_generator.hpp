#ifndef CUT_GENERATOR_HPP
#define CUT_GENERATOR_HPP


#include <linear_program.hpp>

class CutGenerator {
public:
	enum CutStatus {
		valid,
		maybe_recalc,
		recalc
	};

	/**
	 * Prüft, ob die angegebene Lösung gültig ist und fügt, falls sie nicht gültig ist, eine oder mehrere Ungleichungen
	 * zum LP hinzu, die von der Lösung nicht erfüllt werden
	 * @param lp Das LP, zu dem die neuen Ungleichungen hinzugefügt werden sollen
	 * @param solution Die zu betrachtende LP-Lösung
	 * @return valid, falls keine Ungleichungen hinzugefügt wurden
	 * maybe_recalc, falls Ungleichungen hinzugefügt wurden, aber eine Neuberechnung der Lösung nicht unbedingt
	 * notwendig ist, z.B. falls die Ungleichungen nur fraktionale Lösungen entfernen, aber nach ganzzahligen gesucht
	 * wird
	 * recalc, falls Ungleichungen hinzugefügt wurden und die Lösung neu berechnet werden muss
	 */
	virtual CutStatus validate(LinearProgram& lp, const std::vector<double>& solution, CutStatus currentStatus) = 0;
};


#endif