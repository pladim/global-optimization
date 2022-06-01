#ifndef TESTER_H
#define TESTER_H

#include <vector>

#include "SolverFactory.h"
#include "Functions.h"

class Tester {
	shared_ptr<DivideByThree> _solver;
	std::string _name_of_method;
	Parameters _parameters;
	std::vector<uint> _measurements;
	std::vector<double> _distances;
private:
	void solve_1();
	void solve_2();
	void solve_3();
	void solve_4();
	void solve_5();
	void solve_6();
	void solve_gen(const int& task);
public:
	Tester();
	void start_testing();
	std::vector<uint> get_measurements();
	void write_measurements_to_file();
	void solve_task(const int& task);
	std::vector<double> get_distances();
	void write_distances_to_file();
	~Tester();
};

#endif // TESTER_H
