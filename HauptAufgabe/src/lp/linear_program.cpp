#include <cassert>
#include <stdexcept>
#include <linear_program.hpp>
#include <cut_generator.hpp>
#include <cmath>

LinearProgram::LinearProgram(std::string name, Goal opt) : problem(QScreate_prob(name.c_str(), opt)) {
	assert(problem!=nullptr);
}

LinearProgram::~LinearProgram() {
	QSfree_prob(problem);
}

/**
 * Fügt eine neue Variable zum LP hinzu
 * @param objCoeff Der Koeffizient der Variablen in der Zielfunktion
 * @param lower Untere Schranke für die Variable
 * @param upper Obere Schranke für die variable
 */
void LinearProgram::addVariable(double objCoeff, double lower, double upper) {
	int result = QSnew_col(problem, objCoeff, lower, upper, nullptr);
	if (result!=0) {
		throw std::runtime_error("Could not add variable to LP, return value was "+std::to_string(result));
	}
}

/**
 * Fügt eine neue (Un)Gleichung zum LP hinzu
 * @param indices Die Indizes der Koeffizienten, die nicht 0 sind
 * @param coeffs Die Werte dieser Koeffizienten
 * @param rhs Die rechte Seite der Ungleichung
 * @param sense Um welche Art von (Un)Gleichung es sich handelt (kleiner gleich, größer gleich, gleich)
 */
void LinearProgram::addConstraint(const std::vector<int>& indices, const std::vector<double>& coeffs, double rhs,
								  LinearProgram::CompType sense) {
	assert(indices.size()==coeffs.size());
	int result = QSadd_row(problem, static_cast<int>(indices.size()),
						   const_cast<int*>(indices.data()), const_cast<double*>(coeffs.data()),
						   rhs, sense, nullptr);
	if (result!=0) {
		throw std::runtime_error("Could not add constraint to LP, return value was "+std::to_string(result));
	}
}

/**
 * Entfernt die Constraints mit den angegebenen Indizes. Die Indizes der verbleibenden Constraints "rücken auf",
 * d.h. nach dem Entfernen von Constraint 0, 2 und 3 wird Constraint 1 zu Constraint 0, 4 zu 1, 5 zu 2, 6 zu 3, etc.
 * @param indices die Indizes der zu entfernenden Constraints
 */
void LinearProgram::removeConstraints(const std::vector<int>& indices) {
	int status = QSdelete_rows(problem, static_cast<int>(indices.size()), const_cast<int*>(indices.data()));
	if (status!=0) {
		throw std::runtime_error("Error while deleting constraints: "+std::to_string(status));
	}
}

/**
 * @return eine optimale Lösung des LP
 */
LinearProgram::Solution LinearProgram::solve() {
	Solution sol(static_cast<size_t>(getVariableCount()));
	solve(sol);
	return sol;
}

/**
 * Löst das LP. Falls das LP erfolgreich gelöst wurde, wird die gefundene Lösung in out ausgegeben.
 * Falls das LP unzulässig ist, wird out dementsprechend gesetzt. In allen anderen Fällen wird ein
 * Fehler geworfen.
 * @param out Eine Solution-Objekt der korrekten Größe. Dient als Ausgabe.
 */
void LinearProgram::solve(LinearProgram::Solution& out) {
	int status;
	int result = QSopt_dual(problem, &status);
	if (result!=0) {
		throw std::runtime_error("Could not solve LP, return value was "+std::to_string(result));
	}
	switch (status) {
		case QS_LP_OPTIMAL:
			result = QSget_x_array(problem, out.vector.data());
			if (result!=0) {
				throw std::runtime_error("Failed to copy LP solution: "+std::to_string(status));
			}
			result = QSget_objval(problem, &out.value);
			if (result!=0) {
				throw std::runtime_error("Failed to copy LP objective value: "+std::to_string(status));
			}
			break;
		case QS_LP_INFEASIBLE:
			out.value = NAN;
			break;
		case QS_LP_UNBOUNDED:
			throw std::runtime_error("LP is unbounded");
		case QS_LP_ITER_LIMIT:
			throw std::runtime_error("LP solver reached iteration limit");
		case QS_LP_TIME_LIMIT:
			throw std::runtime_error("LP solver reached time limit");
		case QS_LP_UNSOLVED:
			throw std::runtime_error("LP solver failed to solve the problem");
		default:
			throw std::runtime_error("LP solver gave unknown status: "+std::to_string(status));
	}
}

int LinearProgram::getVariableCount() {
	return QSget_colcount(problem);
}

double LinearProgram::getBound(int var, BoundType bound) {
	double ret;
	int result = QSget_bound(problem, var, bound, &ret);
	if (result!=0) {
		throw std::runtime_error("Failed to get variable bound: "+std::to_string(result));
	}
	return ret;
}

void LinearProgram::setBound(int var, LinearProgram::BoundType type, double value) {
	QSchange_bound(problem, var, type, value);
}

LinearProgram::Goal LinearProgram::getGoal() {
	int g;
	QSget_objsense(problem, &g);
	return static_cast<Goal>(g);
}

int LinearProgram::getConstraintCount() {
	return QSget_rowcount(problem);
}

std::vector<double> LinearProgram::getSlack() {
	std::vector<double> ret(getConstraintCount());
	int status = QSget_slack_array(problem, ret.data());
	if (status!=0) {
		throw std::runtime_error("Could not get slack array: "+std::to_string(status));
	}
	return ret;
}

std::vector<double> LinearProgram::getObjective() {
	std::vector<double> ret(getVariableCount());
	int status = QSget_obj(problem, ret.data());
	if (status!=0) {
		throw std::runtime_error("Could not get objective coefficients: "+std::to_string(status));
	}
	return ret;
}

LinearProgram::Solution::Solution(size_t varCount) : vector(varCount) {}

double LinearProgram::Solution::operator[](size_t index) {
	return vector[index];
}

double LinearProgram::Solution::getValue() {
	return value;
}

bool LinearProgram::Solution::isValid() {
	return !std::isnan(value);
}

const std::vector<double>& LinearProgram::Solution::getVector() {
	return vector;
}
