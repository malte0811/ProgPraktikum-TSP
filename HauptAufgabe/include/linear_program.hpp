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

	LinearProgram(std::string name, Goal opt);

	LinearProgram(const LinearProgram& other) = delete;

	~LinearProgram();

	void addVariable(double objCoeff, double lower, double upper);

	void addConstraint(std::vector<int> indices, std::vector<double> coeffs, double rhs,
					   LinearProgram::CompType sense);

	std::vector<double> solve();

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