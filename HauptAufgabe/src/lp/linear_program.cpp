#include <cassert>
#include <stdexcept>
#include <linear_program.hpp>
#include <cmath>
#include <ilcplex/cplex.h>
#include <tsp_utils.hpp>
#include <lemon/tolerance.h>
#include <stddef.h>
#include <string>
#include <utility>
#include <vector>
#include <ilcplex/cpxconst.h>

variable_id LinearProgram::invalid_variable = -1;

LinearProgram::LinearProgram(const SharedCplexEnv& env, const std::string& name, Goal opt) : env(env) {
	int status;
	problem = CPXcreateprob(env.get(), &status, name.c_str());
	if (status != 0) {
		throw std::runtime_error("Could not create CPLEX problem: " + getErrorMessage(status, env.get()));
	}
	CPXchgobjsen(env.get(), problem, opt);
}

LinearProgram::~LinearProgram() {
	CPXfreeprob(env.get(), &problem);
}

/**
 * Fügt neue Variablen zum LP hinzu
 * @param objCoeff Die Koeffizienten der Variablen in der Zielfunktion
 * @param lower Untere Schranken für die Variablen
 * @param upper Obere Schranken für die variablen
 */
void LinearProgram::addVariables(const std::vector<double>& objCoeff, const std::vector<double>& lower,
								 const std::vector<double>& upper) {
	assert(objCoeff.size() == lower.size() && objCoeff.size() == upper.size());
	int result = CPXnewcols(env.get(), problem, static_cast<int>(objCoeff.size()), objCoeff.data(), lower.data(),
							upper.data(),
							nullptr, nullptr);
	if (result != 0) {
		throw std::runtime_error("Could not add variables to LP, return value was " +
								 getErrorMessage(result, env.get()));
	}
	varCount += objCoeff.size();
}

/**
 * Fügt die angegebene Constraint zum LP hinzu
 */
void LinearProgram::addConstraint(const Constraint& constr) {
	std::vector<Constraint> temp{constr};
	addConstraints(temp.begin(), temp.end());
}

/**
 * Entfernt die Constraints, die im Vector den Wert 1 haben und setzt den Wert an Index i entweder auf die neue ID der
 * Constraint i oder auf -1, falls i entfernt wurde
 */
void LinearProgram::removeSetConstraints(std::vector<int>& shouldDelete) {
	assert(constraints.size() == shouldDelete.size());
	int status = CPXdelsetrows(env.get(), problem, shouldDelete.data());
	if (status != 0) {
		throw std::runtime_error("Error while deleting constraints: " + getErrorMessage(status, env.get()));
	}
	std::vector<Constraint> newConstrs;
	for (size_t i = 0; i < constraints.size(); ++i) {
		if (shouldDelete[i] >= 0) {
			newConstrs.push_back(constraints[i]);
		}
	}
	constraints = std::move(newConstrs);
}

/**
 * Entfernt die Variablen mit Wert 1 im vector und setzt den Wert an Index i entweder auf die neue ID der Variablen i
 * oder auf -1, falls i entfernt wurde
 */
void LinearProgram::removeSetVariables(std::vector<int>& shouldDelete) {
	int status = CPXdelsetcols(env.get(), problem, shouldDelete.data());
	for (Constraint& c:constraints) {
		c.deleteVariables(shouldDelete);
	}
	if (status != 0) {
		throw std::runtime_error("Error while deleting variables: " + getErrorMessage(status, env.get()));
	}
}

/**
 * Löst das LP. Falls das LP erfolgreich gelöst wurde, wird die gefundene Lösung in out ausgegeben.
 * Falls das LP unzulässig ist, wird out dementsprechend gesetzt. In allen anderen Fällen wird ein
 * Fehler geworfen.
 * @param out Eine Solution-Objekt der korrekten Größe. Dient als Ausgabe.
 */
void LinearProgram::solveDual(LinearProgram::Solution& out) {
	int result = CPXdualopt(env.get(), problem);
	if (result != 0) {
		throw std::runtime_error("Could not solve LP, return value was " + getErrorMessage(result, env.get()));
	}
	int status = CPXgetstat(env.get(), problem);
	switch (status) {
		case CPX_STAT_OPTIMAL: {
			out.valid = true;
			out.slack.resize(static_cast<size_t>(getConstraintCount()));
			out.shadowCosts.resize(static_cast<size_t>(getConstraintCount()));
			result = CPXsolution(env.get(), problem, &status, &out.value, out.vector.data(), out.shadowCosts.data(),
								 out.slack.data(), out.reduced.data());
			if (result != 0) {
				throw std::runtime_error("Failed to copy LP solution: " + getErrorMessage(result, env.get()));
			}
		}
			break;
		case CPX_STAT_INFEASIBLE:
		case CPX_STAT_INForUNBD:
			out.valid = false;
			break;
		case CPX_STAT_UNBOUNDED:
			throw std::runtime_error("LP is unbounded");
		default:
			throw std::runtime_error("LP solver gave error: " + getErrorMessage(status, env.get()));
	}
}

variable_id LinearProgram::getVariableCount() {
	return varCount;
}

double LinearProgram::getBound(variable_id var, BoundType bound) {
	double ret;
	int result;
	if (bound == lower) {
		result = CPXgetlb(env.get(), problem, &ret, var, var);
	} else {
		result = CPXgetub(env.get(), problem, &ret, var, var);
	}
	if (result != 0) {
		throw std::runtime_error("Failed to get variable bound: " + getErrorMessage(result, env.get()));
	}
	return ret;
}

void LinearProgram::setBound(variable_id var, LinearProgram::BoundType type, double value) {
	char typeChar = static_cast<char>(type);
	CPXchgbds(env.get(), problem, 1, &var, &typeChar, &value);
}

LinearProgram::Goal LinearProgram::getGoal() {
	return static_cast<Goal>(CPXgetobjsen(env.get(), problem));
}

int LinearProgram::getConstraintCount() {
	return constraints.size();
}

std::vector<double> LinearProgram::getObjective() {
	std::vector<double> ret(static_cast<size_t>(varCount));
	int status = CPXgetobj(env.get(), problem, ret.data(), 0, varCount - 1);
	if (status != 0) {
		throw std::runtime_error("Could not get objective coefficients: " + getErrorMessage(status, env.get()));
	}
	return ret;
}

LinearProgram::Constraint LinearProgram::getConstraint(int index) const {
	return constraints[index];
}

SharedCplexEnv LinearProgram::openCPLEX() {
	int status;
	SharedCplexEnv ret(CPXopenCPLEX(&status), [](CPXENVptr& env) {
		CPXcloseCPLEX(&env);
	});
	if (status != 0) {
		throw std::runtime_error("Failed to open CPLEX environment: " +
								 LinearProgram::getErrorMessage(status, nullptr));
	}
	return ret;
}

std::string LinearProgram::getErrorMessage(int error, CPXENVptr env) {
	char buffer[4096];
	if (CPXgeterrorstring(env, error, buffer) != nullptr) {
		return buffer;
	}
	return "Unkown error: " + std::to_string(error);
}

LinearProgram::Solution::Solution(size_t varCount, size_t constraintCount)
		: vector(varCount), slack(constraintCount), reduced(varCount), shadowCosts(constraintCount), value(0),
		  valid(true) {}

double LinearProgram::Solution::operator[](size_t index) const {
	return vector[index];
}

double LinearProgram::Solution::getValue() const {
	return value;
}

bool LinearProgram::Solution::isValid() const {
	return valid;
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

const std::vector<double>& LinearProgram::Solution::getShadowCosts() const {
	return shadowCosts;
}

void LinearProgram::Solution::removeVariables(const std::vector<variable_id>& toRemove) {
	tsp_util::eraseEntries(vector, toRemove);
	tsp_util::eraseEntries(reduced, toRemove);
}

LinearProgram::Constraint::Constraint(const std::vector<int>& indices, const std::vector<double>& coeffs,
									  LinearProgram::CompType cmp, double rhs) :
		indices(indices), coeffs(coeffs), comp(cmp), rhs(rhs) {
	assert(indices.size() == coeffs.size());
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
	for (size_t i = 0; i < indices.size(); ++i) {
		ret += variables[indices[i]] * coeffs[i];
	}
	return ret;
}

LinearProgram::Constraint::Constraint() : Constraint({}, {}, less_eq, -1) {}

bool LinearProgram::Constraint::isValid() const {
	return !indices.empty();
}

void LinearProgram::Constraint::deleteVariables(const std::vector<variable_id>& removalMap) {
	std::vector<variable_id> newNZ;
	std::vector<double> newCoeffs;
	for (size_t i = 0; i < coeffs.size(); ++i) {
		if (removalMap[indices[i]] >= 0) {
			newNZ.push_back(removalMap[indices[i]]);
			newCoeffs.push_back(coeffs[i]);
		}
	}
	coeffs = newCoeffs;
	indices = newNZ;
	assert(!indices.empty());
	assert(indices.size() == coeffs.size());
}

bool LinearProgram::Constraint::isViolated(const std::vector<double>& vars, lemon::Tolerance<double> tolerance) const {
	return !isValidLHS(evalLHS(vars), tolerance);
}

LinearProgram::Constraint LinearProgram::Constraint::fromDense(const std::vector<variable_id>& usedVars,
															   const std::vector<double>& coeffs,
															   LinearProgram::CompType sense, double rhs) {
	std::vector<double> usedCoeffs;
	usedCoeffs.reserve(usedVars.size());
	for (variable_id var:usedVars) {
		usedCoeffs.push_back(coeffs[var]);
	}
	return Constraint(usedVars, usedCoeffs, sense, rhs);
}