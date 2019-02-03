#ifndef LINEAR_PROGRAM_HPP
#define LINEAR_PROGRAM_HPP

#include <vector>
#include <qsopt.h>
#include <string>

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

		double operator[](size_t index);

		double getValue();

		bool isValid();

		const std::vector<double>& getVector();

	private:
		friend LinearProgram;
		std::vector<double> vector;
		double value;
	};

	LinearProgram(std::string name, Goal opt);

	LinearProgram(const LinearProgram& other) = delete;

	~LinearProgram();

	void addVariable(double objCoeff, double lower, double upper);

	void addConstraint(std::vector<int> indices, std::vector<double> coeffs, double rhs,
					   LinearProgram::CompType sense);

	Solution solve();

	void solve(Solution& out);

	static constexpr double pos_infinite = QS_MAXDOUBLE;
	static constexpr double neg_infinite = -QS_MAXDOUBLE;

	int getVariableCount();

	double getBound(int var, BoundType bound);

	void setBound(int var, BoundType type, double value);

	Goal getGoal();

private:
	QSprob problem;
};


#endif