#ifndef LINEAR_PROGRAM_HPP
#define LINEAR_PROGRAM_HPP

#include <vector>
#include <string>
#include <ilcplex/cplex.h>

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
		explicit Solution(size_t varCount);

		double operator[](size_t index) const;

		double getValue() const;

		bool isValid() const;

		const std::vector<double>& getVector() const;

		const std::vector<double>& getSlack() const;

		const std::vector<double>& getReducedCosts() const;

	private:
		friend LinearProgram;
		std::vector<double> vector;
		std::vector<double> slack;
		std::vector<double> reduced;
		double value;
	};

	LinearProgram(CPXENVptr& env, std::string name, Goal opt);

	LinearProgram(const LinearProgram& other) = delete;

	~LinearProgram();

	void addVariables(const std::vector<double>& objCoeff, const std::vector<double>& lower,
					  const std::vector<double>& upper);

	void addConstraint(const std::vector<int>& indices, const std::vector<double>& coeffs, double rhs,
					   LinearProgram::CompType sense);

	void removeSetConstraints(std::vector<int>& indices);

	Solution solve();

	void solve(Solution& out);

	int getVariableCount();

	int getConstraintCount();

	double getBound(int var, BoundType bound);

	void setBound(int var, BoundType type, double value);

	Goal getGoal();

	std::vector<double> getObjective();

private:
	CPXENVptr& env;
	CPXLPptr problem;
};


#endif