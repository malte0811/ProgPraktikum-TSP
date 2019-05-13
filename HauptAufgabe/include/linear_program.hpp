#ifndef LINEAR_PROGRAM_HPP
#define LINEAR_PROGRAM_HPP

#include <vector>
#include <string>
#include <ilcplex/cplex.h>
#include <lemon/tolerance.h>
#include <cassert>
#include <stdexcept>
#include <set>
#include <memory>
#include <type_traits>

using variable_id = int;

using SharedCplexEnv = std::shared_ptr<std::remove_pointer<CPXENVptr>::type>;

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
		bool valid;
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

	LinearProgram(const SharedCplexEnv& env, const std::string& name, Goal opt);

	LinearProgram(const LinearProgram& other) = delete;

	~LinearProgram();

	void addVariables(const std::vector<double>& objCoeff, const std::vector<double>& lower,
					  const std::vector<double>& upper);

	void addConstraint(const Constraint& constr);

	template<typename It>
	void addConstraints(const It& begin, const It& end);

	Constraint getConstraint(int index) const;

	void removeSetConstraints(std::vector<int>& shouldDelete);

	void removeSetVariables(std::vector<int>& shouldDelete);

	void solveDual(Solution& out);

	variable_id getVariableCount();

	int getConstraintCount();

	double getBound(variable_id var, BoundType bound);

	void setBound(variable_id var, BoundType type, double value);

	Goal getGoal();

	std::vector<double> getObjective();

	static variable_id invalid_variable;

	static SharedCplexEnv openCPLEX();
private:
	static std::string getErrorMessage(int error, CPXENVptr env);

	SharedCplexEnv env;
	CPXLPptr problem;
	std::vector<Constraint> constraints;
	variable_id varCount = 0;
};

template<typename It>
void LinearProgram::addConstraints(const It& begin, const It& end) {
	std::vector<double> rhs;
	std::vector<double> coeffs;
	std::vector<int> indices;
	std::vector<int> constrStarts;
	std::vector<char> sense;
	size_t addCount = std::distance(begin, end);
	constrStarts.reserve(addCount);
	sense.reserve(addCount);
	rhs.reserve(addCount);
	constraints.reserve(constraints.size() + addCount);
	for (It i = begin; i != end; ++i) {
		const auto& c = static_cast<const Constraint&>(*i);
		assert(c.isValid() || getVariableCount() == 0);
		constrStarts.push_back(indices.size());
		rhs.push_back(c.getRHS());
		sense.push_back(c.getSense());
		indices.insert(indices.end(), c.getNonzeroes().begin(), c.getNonzeroes().end());
		coeffs.insert(coeffs.end(), c.getCoeffs().begin(), c.getCoeffs().end());
		constraints.push_back(c);
	}
	int result = CPXaddrows(env.get(), problem, 0, addCount, indices.size(), rhs.data(),
							sense.data(), constrStarts.data(), indices.data(), coeffs.data(), nullptr, nullptr);
	if (result != 0) {
		throw std::runtime_error("Could not add constraint to LP, return value was " +
								 getErrorMessage(result, env.get()));
	}
}

#endif