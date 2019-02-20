#ifndef LINEAR_PROGRAM_HPP
#define LINEAR_PROGRAM_HPP

#include <vector>
#include <qsopt.h>
#include <string>

using variable_id = int;

class LinearProgram {
public:
	enum CompType {
		less_eq = 'L',
		equal = 'E',
		greater_eq = 'G'
	};
	enum Goal {
		minimize = QS_MIN,
		maximize = QS_MAX
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

	LinearProgram(std::string name, Goal opt);

	LinearProgram(const LinearProgram& other) = delete;

	~LinearProgram();

	void addVariable(double objCoeff, double lower, double upper);

	void addConstraint(const std::vector<int>& indices, const std::vector<double>& coeffs, double rhs,
					   LinearProgram::CompType sense);

	void removeConstraints(const std::vector<int>& indices);

	Solution solve();

	void solve(Solution& out);

	static constexpr double pos_infinite = QS_MAXDOUBLE;
	static constexpr double neg_infinite = -QS_MAXDOUBLE;

	int getVariableCount();

	int getConstraintCount();

	double getBound(int var, BoundType bound);

	void setBound(int var, BoundType type, double value);

	Goal getGoal();

	std::vector<double> getObjective();

private:
	QSprob problem;
};


#endif