#ifndef LINEAR_PROGRAM_HPP
#define LINEAR_PROGRAM_HPP

#include <vector>
#include <string>
#include <ilcplex/cplex.h>
#include <lemon/tolerance.h>
#include <cassert>
#include <stdexcept>
#include <set>

using variable_id = int;

class LinearProgram {
public:
	enum CompType {
		less_eq = 'L',
		equal = 'E',
		greater_eq = 'G'
	};
	enum Goal {
		minimize = CPX_MIN,
		maximize = CPX_MAX
	};
	enum BoundType {
		lower = 'L',
		upper = 'U'
	};

	class Solution {
	public:
		Solution(size_t varCount, size_t constraintCount);

		double operator[](size_t index) const;

		double getValue() const;

		bool isValid() const;

		const std::vector<double>& getVector() const;

		const std::vector<double>& getSlack() const;

		const std::vector<double>& getReducedCosts() const;

		const std::vector<double>& getShadowCosts() const;

		void removeVariables(const std::vector<variable_id>& toRemove);

	private:
		friend LinearProgram;
		std::vector<double> vector;
		std::vector<double> slack;
		std::vector<double> reduced;
		std::vector<double> shadowCosts;
		double value;
	};

	class Constraint {
	public:
		Constraint(const std::vector<int>& indices, const std::vector<double>& coeffs, CompType cmp, double rhs);

		Constraint();

		const std::vector<int>& getNonzeroes() const;

		const std::vector<double>& getCoeffs() const;

		CompType getSense() const;

		double getRHS() const;

		bool isValidLHS(double lhs, lemon::Tolerance<double> tolerance) const;

		double evalLHS(const std::vector<double>& variables) const;

		bool isValid() const;

		bool isViolated(const std::vector<double>& vars, lemon::Tolerance<double> tolerance) const;

		void deleteVariables(const std::vector<variable_id>& removalMap);

		static Constraint fromDense(const std::vector<variable_id>& usedVars, const std::vector<double>& coeffs,
									CompType sense,
									double rhs);

	private:
		std::vector<int> indices;
		std::vector<double> coeffs;
		CompType comp;
		double rhs;
		friend LinearProgram;
	};

	LinearProgram(CPXENVptr& env, const std::string& name, Goal opt);

	LinearProgram(const LinearProgram& other) = delete;

	~LinearProgram();

	void addVariables(const std::vector<double>& objCoeff, const std::vector<double>& lower,
					  const std::vector<double>& upper);

	void addVariableWithCoeffs(double objCoeff, double lower, double upper, const std::vector<int>& indices,
							   const std::vector<double>& constrCoeffs);

	void addConstraint(const Constraint& constr);

	template<typename T>
	void addConstraints(const T& constrs);

	Constraint getConstraint(int index) const;

	void removeSetConstraints(std::vector<int>& shouldDelete);

	void removeSetVariables(std::vector<int>& shouldDelete);

	Solution solveDual();

	Solution solvePrimal();

	void solve(Solution& out);

	variable_id getVariableCount();

	int getConstraintCount();

	double getBound(variable_id var, BoundType bound);

	void setBound(variable_id var, BoundType type, double value);

	Goal getGoal();

	std::vector<double> getObjective();

	static variable_id invalid_variable;
private:
	void writeSolution(Solution& sol);

	std::string getErrorMessage(int error);

	CPXENVptr& env;
	CPXLPptr problem;
	std::vector<Constraint> constraints;
};

template<typename T>
void LinearProgram::addConstraints(const T& constrs) {
	std::vector<double> rhs;
	std::vector<double> coeffs;
	std::vector<int> indices;
	std::vector<int> constrStarts;
	std::vector<char> sense;
	constrStarts.reserve(constrs.size());
	sense.reserve(constrs.size());
	rhs.reserve(constrs.size());
	for (const Constraint& c:constrs) {
		assert(c.isValid() || getVariableCount() == 0);
		constrStarts.push_back(indices.size());
		rhs.push_back(c.getRHS());
		sense.push_back(c.getSense());
		indices.insert(indices.end(), c.getNonzeroes().begin(), c.getNonzeroes().end());
		coeffs.insert(coeffs.end(), c.getCoeffs().begin(), c.getCoeffs().end());
	}
	int result = CPXaddrows(env, problem, 0, constrs.size(), indices.size(), rhs.data(),
							sense.data(), constrStarts.data(), indices.data(), coeffs.data(), nullptr, nullptr);
	if (result != 0) {
		throw std::runtime_error("Could not add constraint to LP, return value was " + getErrorMessage(result));
	}
	constraints.insert(constraints.end(), constrs.begin(), constrs.end());
}

#endif