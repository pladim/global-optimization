#include <fstream>

#include "Tester.h"
#include "Functions.h"
#include "SolverFactory.h"

Tester::Tester(const std::string& name_of_method,
			   const Parameters& parameters) :
	_name_of_method(name_of_method),
	_parameters(parameters) {
	init_dll();

	GKLSG_Int_Void __GetDim = (GKLSG_Int_Void)GetProcAddress(hDLL, "GetDim");
	GKLSG_Int_Void __GetCountCondition = (GKLSG_Int_Void)GetProcAddress(hDLL, "GetCountCondition");

	uint dim{ static_cast<uint>(__GetDim()) };
	uint cst{ static_cast<uint>(__GetCountCondition()) };

	_parameters._dimension = dim;
	_parameters._constraints = cst;
	_parameters._eps = sqrtl(static_cast<double>(dim)) * 1.0 + 0.5;
}

void Tester::start_testing() {
	if (_parameters._state != Mode::test) return;

	for (int i = 1; i < 101; ++i) {
		generate_task(i);
		std::cout << "\t\t\tTASK #" << i << std::endl;
		CoordinatesValues solution_point{ get_minimum_point() };
		double solution_minimum = global_minimum(solution_point);
		std::cout << "True global minimum is " << solution_minimum << std::endl;
		std::cout << "was reached at the point:" << std::endl;

		for (auto& elem : solution_point)
			std::cout << elem << std::endl;

		Problem testDll(_parameters._dimension, _parameters._constraints, { -1, -1 }, { 1, 1 }, &test_dll);
		testDll.set_true_minimum(solution_point);

		_solver = create_solver(
			_name_of_method, 
			_parameters._dimension, 
			_parameters._constraints, 
			_parameters, 
			testDll
		);

		_solver->solve();
		if (_solver->solved()) _measurements.push_back(_solver->get_gen());

		std::cout << std::endl;
	}

	write_measurements_to_file();
}

std::vector<uint> Tester::get_measurements() {
	return _measurements;
}

void Tester::write_measurements_to_file() {
	if (_parameters._state == Mode::test) {
		std::string _path{ "..\\x64\\Release\\" };
		_path = _path + _name_of_method + "_measurements.txt";
		std::ofstream out;
		out.open(_path);

		if (out.is_open()) {
			for (auto& elem : _measurements)
				out << elem << std::endl;
		}
	}
}

void Tester::solve_1() {
	std::cout << "\t\t\tTASK #1 (NON-GEN)" << std::endl;
	std::cout << "True global minimum is -1.48968" << std::endl;
	std::cout << "was reached at the point:" << std::endl;
	std::cout << "0.94248" << std::endl;
	std::cout << "0.94526" << std::endl;

	uint dim_t1{ 2 };
	uint cst_t1{ 3 };
	Problem testProblem1(dim_t1, cst_t1, { 0.0, -1.0 }, { 4.0, 3.0 }, &task1);

	uint dm = _parameters._dimension;
	uint ct = _parameters._constraints;
	_parameters._dimension = dim_t1;
	_parameters._constraints = cst_t1;
	_parameters._eps = sqrtl(static_cast<double>(dim_t1)) * 1.0 + 0.5;

	_solver = create_solver(_name_of_method, dim_t1, cst_t1, _parameters, testProblem1);

	_solver->solve();
	_solver->write_generated_points();
	_solver->write_generated_intervals();

	_distances = _solver->get_distances();
	write_distances_to_file();

	_parameters._dimension = dm;
	_parameters._constraints = ct;
}

void Tester::solve_2() {
	std::cout << "\t\t\tTASK #2 (NON-GEN)" << std::endl;
	std::cout << "True global minimum is -0.80467" << std::endl;
	std::cout << "was reached at the point:" << std::endl;
	std::cout << "-0.32520" << std::endl;
	std::cout << "0.78197" << std::endl;

	uint dim_t2{ 2 };
	uint cst_t2{ 2 };
	Problem testProblem2(dim_t2, cst_t2, { -2.5, -1.5 }, { 2.5, 1.5 }, &task2);

	uint dm = _parameters._dimension;
	uint ct = _parameters._constraints;
	_parameters._dimension = dim_t2;
	_parameters._constraints = cst_t2;
	_parameters._eps = sqrtl(static_cast<double>(dim_t2)) * 1.0 + 0.5;

	_solver = create_solver(_name_of_method, dim_t2, cst_t2, _parameters, testProblem2);

	_solver->solve();
	_solver->write_generated_points();
	_solver->write_generated_intervals();

	_distances = _solver->get_distances();
	write_distances_to_file();

	_parameters._dimension = dm;
	_parameters._constraints = ct;
}

void Tester::solve_3() {
	std::cout << "\t\t\tTASK #3 (NON-GEN)" << std::endl;
	std::cout << "True global minimum is -0.81911" << std::endl;
	std::cout << "was reached at the point:" << std::endl;
	std::cout << "1.30499" << std::endl;
	std::cout << "2.27249" << std::endl;

	uint dim_t3{ 2 };
	uint cst_t3{ 1 };
	Problem testProblem3(dim_t3, cst_t3, { 0, 0 }, { 2 * _PI, 2 * _PI }, &task3);

	uint dm = _parameters._dimension;
	uint ct = _parameters._constraints;
	_parameters._dimension = dim_t3;
	_parameters._constraints = cst_t3;
	_parameters._eps = sqrtl(static_cast<double>(dim_t3)) * 1.0 + 0.5;

	_solver = create_solver(_name_of_method, dim_t3, cst_t3, _parameters, testProblem3);

	_solver->solve();
	_solver->write_generated_points();
	_solver->write_generated_intervals();

	_distances = _solver->get_distances();
	write_distances_to_file();

	_parameters._dimension = dm;
	_parameters._constraints = ct;
}

void Tester::solve_4() {
	std::cout << "\t\t\tTASK #4 (NON-GEN)" << std::endl;
	std::cout << "True global minimum is -1.97384" << std::endl;
	std::cout << "was reached at the point:" << std::endl;
	std::cout << "0.1" << std::endl;
	std::cout << "0.0" << std::endl;

	uint dim_t4{ 2 };
	uint cst_t4{ 1 };
	Problem testProblem4(dim_t4, cst_t4, { -1.0, -1.0 }, { 1.0, 1.0 }, &task4);

	uint dm = _parameters._dimension;
	uint ct = _parameters._constraints;
	_parameters._dimension = dim_t4;
	_parameters._constraints = cst_t4;
	_parameters._eps = sqrtl(static_cast<double>(dim_t4)) * 1.0 + 0.5;

	_solver = create_solver(_name_of_method, dim_t4, cst_t4, _parameters, testProblem4);

	_solver->solve();
	_solver->write_generated_points();
	_solver->write_generated_intervals();

	_distances = _solver->get_distances();
	write_distances_to_file();

	_parameters._dimension = dm;
	_parameters._constraints = ct;
}

void Tester::solve_5() {
	std::cout << "\t\t\tTASK #5 (NON-GEN)" << std::endl;
	std::cout << "True global minimum is 0.09768" << std::endl;
	std::cout << "was reached at the point:" << std::endl;
	std::cout << "0.9" << std::endl;
	std::cout << "1.0" << std::endl;

	uint dim_t5{ 2 };
	uint cst_t5{ 1 };
	Problem testProblem5(dim_t5, cst_t5, { -10, -10 }, { 10, 10 }, &task5);

	uint dm = _parameters._dimension;
	uint ct = _parameters._constraints;
	_parameters._dimension = dim_t5;
	_parameters._constraints = cst_t5;
	_parameters._eps = sqrtl(static_cast<double>(dim_t5)) * 1.0 + 1.5;

	_solver = create_solver(_name_of_method, dim_t5, cst_t5, _parameters, testProblem5);

	_solver->solve();
	_solver->write_generated_points();
	_solver->write_generated_intervals();

	_distances = _solver->get_distances();
	write_distances_to_file();

	_parameters._dimension = dm;
	_parameters._constraints = ct;
}

void Tester::solve_6() {
	std::cout << "\t\t\tTASK #6 (NON-GEN)" << std::endl;
	std::cout << "True global minimum is -0.81911" << std::endl;
	std::cout << "was reached at the point:" << std::endl;
	std::cout << "1.30499" << std::endl;
	std::cout << "2.27249" << std::endl;

	uint dim_t6{ 2 };
	uint cst_t6{ 1 };
	Problem testProblem6(dim_t6, cst_t6, { 0, 0 }, { 2 * _PI, 2 * _PI }, &task6);

	uint dm = _parameters._dimension;
	uint ct = _parameters._constraints;
	_parameters._dimension = dim_t6;
	_parameters._constraints = cst_t6;
	_parameters._eps = sqrtl(static_cast<double>(dim_t6)) * 1.0 + 0.5;

	_solver = create_solver(_name_of_method, dim_t6, cst_t6, _parameters, testProblem6);

	_solver->solve();
	_solver->write_generated_points();
	_solver->write_generated_intervals();

	_distances = _solver->get_distances();
	write_distances_to_file();

	_parameters._dimension = dm;
	_parameters._constraints = ct;
}

void Tester::solve_gen(const int& task) {
	generate_task(task);
	std::cout << "\t\t\tTASK #" << task << std::endl;
	CoordinatesValues solution_point{ get_minimum_point() };
	double solution_minimum = global_minimum(solution_point);
	std::cout << "True global minimum is " << solution_minimum << std::endl;
	std::cout << "was reached at the point:" << std::endl;

	for (auto& elem : solution_point)
		std::cout << elem << std::endl;

	GKLSG_Int_Void __GetDim = (GKLSG_Int_Void)GetProcAddress(hDLL, "GetDim");
	GKLSG_Int_Void __GetCountCondition = (GKLSG_Int_Void)GetProcAddress(hDLL, "GetCountCondition");

	uint dim{ static_cast<uint>(__GetDim()) };
	uint cst{ static_cast<uint>(__GetCountCondition()) };

	uint dm = _parameters._dimension;
	uint ct = _parameters._constraints;
	_parameters._dimension = dim;
	_parameters._constraints = cst;
	_parameters._eps = sqrtl(static_cast<double>(dim)) * 1.0 + 0.5;

	Problem testDll(dim, cst, { -1, -1 }, { 1, 1 }, &test_dll);
	_solver = create_solver(_name_of_method, dim, cst, _parameters, testDll);

	_solver->solve();
	_solver->write_generated_points();
	_solver->write_generated_intervals();

	_distances = _solver->get_distances();
	write_distances_to_file();

	_parameters._dimension = dm;
	_parameters._constraints = ct;
}

void Tester::solve_task(const int& task) {
	if ((task > 0) && (task < 101)) {
		solve_gen(task);
	}
	else {
		switch (task % 100) {
		case 1:
			solve_1();
			break;
		case 2:
			solve_2();
			break;
		case 3:
			solve_3();
			break;
		case 4:
			solve_4();
			break;
		case 5:
			solve_5();
			break;
		case 6:
			solve_6();
			break;
		default:
			break;
		}
	}
}

std::vector<double> Tester::get_distances() {
	if (_parameters._state == Mode::stop_by_precision)
		return _distances;
}

void Tester::write_distances_to_file() {
	std::string _path{ "..\\x64\\Release\\" };
	_path = _path + _name_of_method + "_distances.txt";
	std::ofstream out;
	out.open(_path);

	if (out.is_open()) {
		for (auto& elem : _distances)
			out << elem << std::endl;
	}
}

Tester::~Tester() {
	deinit_dll();
}