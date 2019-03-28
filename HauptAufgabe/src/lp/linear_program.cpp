#include <cassert>
#include <stdexcept>
#include <linear_program.hpp>
#include <cut_generator.hpp>
#include <cmath>
#include <ilcplex/cplex.h>

variable_id LinearProgram::invalid_variable = -1;

LinearProgram::LinearProgram(CPXENVptr& env, std::string name, Goal opt) : env(env) {
	int status;
	problem = CPXcreateprob(env, &status, name.c_str());
	if (status!=0) {
		throw std::runtime_error("Could not create CPLEX problem: "+std::to_string(status));
	}
	CPXchgobjsen(env, problem, opt);
}

LinearProgram::~LinearProgram() {
	CPXfreeprob(env, &problem);
}

/**
 * Fügt eine neue Variable zum LP hinzu
 * @param objCoeff Der Koeffizient der Variablen in der Zielfunktion
 * @param lower Untere Schranke für die Variable
 * @param upper Obere Schranke für die variable
 */
void LinearProgram::addVariables(const std::vector<double>& objCoeff, const std::vector<double>& lower,
								 const std::vector<double>& upper) {
	assert(objCoeff.size()==lower.size() && objCoeff.size()==upper.size());
	int result = CPXnewcols(env, problem, static_cast<int>(objCoeff.size()), objCoeff.data(), lower.data(),
							upper.data(),
							nullptr, nullptr);
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
void LinearProgram::addConstraint(const Constraint& constr) {
	addConstraints({constr});
}

void LinearProgram::addConstraints(const std::vector<Constraint>& constrs) {
	std::vector<double> rhs;
	std::vector<double> coeffs;
	std::vector<int> indices;
	std::vector<int> constrStarts;
	std::vector<char> sense;
	constrStarts.reserve(constrs.size());
	sense.reserve(constrs.size());
	rhs.reserve(constrs.size());
	for (const Constraint& c:constrs) {
		assert(c.isValid());
		constrStarts.push_back(indices.size());
		rhs.push_back(c.getRHS());
		sense.push_back(c.getSense());
		indices.insert(indices.end(), c.getNonzeroes().begin(), c.getNonzeroes().end());
		coeffs.insert(coeffs.end(), c.getCoeffs().begin(), c.getCoeffs().end());
	}
	int result = CPXaddrows(env, problem, 0, constrs.size(), indices.size(), rhs.data(),
							sense.data(), constrStarts.data(), indices.data(), coeffs.data(), nullptr, nullptr);
	if (result!=0) {
		throw std::runtime_error("Could not add constraint to LP, return value was "+std::to_string(result));
	}
	constraints.insert(constraints.end(), constrs.begin(), constrs.end());
}

/**
 * Entfernt die Constraints, die im Vector den Wert 1 haben
 */
void LinearProgram::removeSetConstraints(std::vector<int>& shouldDelete) {
	assert(constraints.size()==shouldDelete.size());
	int status = CPXdelsetrows(env, problem, shouldDelete.data());
	if (status!=0) {
		throw std::runtime_error("Error while deleting constraints: "+std::to_string(status));
	}
	std::vector<Constraint> newConstrs;
	for (size_t i = 0; i<constraints.size(); ++i) {
		if (shouldDelete[i]>=0) {
			newConstrs.push_back(constraints[i]);
		}
	}
	constraints = std::move(newConstrs);
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
	int result = CPXdualopt(env, problem);
	if (result!=0) {
		throw std::runtime_error("Could not solve LP, return value was "+std::to_string(result));
	}
	int status = CPXgetstat(env, problem);
	switch (status) {
		case CPX_STAT_OPTIMAL:
			out.slack.resize(static_cast<size_t>(getConstraintCount()));
			//TODO ist nullptr ok?
			result = CPXsolution(env, problem, &status, &out.value, out.vector.data(), nullptr, out.slack.data(),
								 out.reduced.data());
			if (result!=0) {
				throw std::runtime_error("Failed to copy LP solution: "+std::to_string(status));
			}
			break;
		case CPX_STAT_INFEASIBLE:
		case CPX_STAT_INForUNBD:
			//TODO darf ich davon ausgehen, dass NAN existiert/definiert ist?
			out.value = NAN;
			break;
		case CPX_STAT_UNBOUNDED:
			throw std::runtime_error("LP is unbounded");
		default:
			char errStr[510];
			CPXgetstatstring(env, status, errStr);
			throw std::runtime_error("LP solver gave error: "+std::string(errStr));
	}
}

variable_id LinearProgram::getVariableCount() {
	return CPXgetnumcols(env, problem);
}

double LinearProgram::getBound(variable_id var, BoundType bound) {
	double ret;
	int result;
	if (bound==lower) {
		result = CPXgetlb(env, problem, &ret, var, var);
	} else {
		result = CPXgetub(env, problem, &ret, var, var);
	}
	if (result!=0) {
		throw std::runtime_error("Failed to get variable bound: "+std::to_string(result));
	}
	return ret;
}

void LinearProgram::setBound(variable_id var, LinearProgram::BoundType type, double value) {
	char typeChar = static_cast<char>(type);
	CPXchgbds(env, problem, 1, &var, &typeChar, &value);
}

LinearProgram::Goal LinearProgram::getGoal() {
	return static_cast<Goal>(CPXgetobjsen(env, problem));
}

int LinearProgram::getConstraintCount() {
	return constraints.size();
}

std::vector<double> LinearProgram::getObjective() {
	variable_id varCount = getVariableCount();
	std::vector<double> ret(static_cast<size_t>(varCount));
	int status = CPXgetobj(env, problem, ret.data(), 0, varCount-1);
	if (status!=0) {
		throw std::runtime_error("Could not get objective coefficients: "+std::to_string(status));
	}
	return ret;
}

LinearProgram::Constraint LinearProgram::getConstraint(int index) const {
	return constraints[index];
}

LinearProgram::Solution::Solution(size_t varCount) : vector(varCount), slack(0), reduced(varCount) {}

double LinearProgram::Solution::operator[](size_t index) const {
	return vector[index];
}

double LinearProgram::Solution::getValue() const {
	return value;
}

bool LinearProgram::Solution::isValid() const {
	return !std::isnan(value);
}

const std::vector<double>& LinearProgram::Solution::getVector() const {
	return vector;
}

const std::vector<double>& LinearProgram::Solution::getReducedCosts() const {
	return reduced;
}

const std::vector<double>& LinearProgram::Solution::getSlack() const {
	return slack;
}

LinearProgram::Constraint::Constraint(const std::vector<int>& indices, const std::vector<double>& coeffs,
									  LinearProgram::CompType cmp, double rhs) :
		indices(indices), coeffs(coeffs), comp(cmp), rhs(rhs) {
	assert(indices.size()==coeffs.size());
}

const std::vector<int>& LinearProgram::Constraint::getNonzeroes() const {
	return indices;
}

const std::vector<double>& LinearProgram::Constraint::getCoeffs() const {
	return coeffs;
}

LinearProgram::CompType LinearProgram::Constraint::getSense() const {
	return comp;
}

double LinearProgram::Constraint::getRHS() const {
	return rhs;
}

bool LinearProgram::Constraint::isValidLHS(double lhs, lemon::Tolerance<double> tolerance) const {
	switch (comp) {
		case less_eq:
			return !tolerance.less(rhs, lhs);
		case equal:
			return !tolerance.different(lhs, rhs);
		case greater_eq:
			return !tolerance.less(lhs, rhs);
	}
	return false;
}

double LinearProgram::Constraint::evalLHS(const std::vector<double>& variables) const {
	double ret = 0;
	for (size_t i = 0; i<indices.size(); ++i) {
		ret += variables[indices[i]]*coeffs[i];
	}
	return ret;
}

LinearProgram::Constraint::Constraint() : Constraint({}, {}, less_eq, -1) {}

bool LinearProgram::Constraint::isValid() const {
	return !indices.empty();
}
