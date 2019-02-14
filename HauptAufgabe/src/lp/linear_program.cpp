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

void LinearProgram::addVariable(double objCoeff, double lower, double upper) {
	int emptyInt[1];
	double emptyDouble[1];
	int result = QSadd_col(problem, 0, emptyInt, emptyDouble, objCoeff, lower, upper, nullptr);
	if (result!=0) {
		throw std::runtime_error("Could not add variable to LP, return value was "+std::to_string(result));
	}
}

void LinearProgram::addConstraint(std::vector<int> indices, std::vector<double> coeffs, double rhs,
								  LinearProgram::CompType sense) {
	assert(indices.size()==coeffs.size());
	int result = QSadd_row(problem, static_cast<int>(indices.size()), indices.data(), coeffs.data(), rhs, sense,
						   nullptr);
	if (result!=0) {
		throw std::runtime_error("Could not add constraint to LP, return value was "+std::to_string(result));
	}
}

void LinearProgram::removeConstraints(std::vector<int>& indices) {
	int status = QSdelete_rows(problem, static_cast<int>(indices.size()), indices.data());
	if (status!=0) {
		throw std::runtime_error("Error while deleting constraints: "+std::to_string(status));
	}
}

LinearProgram::Solution LinearProgram::solve() {
	Solution sol(static_cast<size_t>(getVariableCount()));
	solve(sol);
	return sol;
}

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
