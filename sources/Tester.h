#ifndef TESTER_H
#define TESTER_H

#include "SolverFactory.h"

class Tester {
	size_t _limit_iter;
	std::string _method;
	Parameters _param;
	shared_ptr<DivideByThree> _solver;
	shared_ptr<Problem> _problem;
	std::vector<double> _solution_point;
	std::vector<double> _numerical_point;
	std::vector<double> _measure_solved;
	double _Delta;
private:
	void run_test(const size_t& max_iter);
public:
	Tester(const std::string& method, Parameters& param, const size_t& limit_iter, const double& Delta);
	void run();
	void write_results();
	~Tester();
};

#endif TESTER_H
