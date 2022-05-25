#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <unordered_map>

#include "synonymous_types.h"
#include "Functions.h"
#include "Problem.h"
#include "Parameters.h"
#include "DivideByThree.h"
#include "SolverFactory.h"

int main() {
	init_dll();

	uint dim_t1{ 2 };
	uint cst_t1{ 3 };
	uint dim_t2{ 2 };
	uint cst_t2{ 2 };
	uint dim_t3{ 2 };
	uint cst_t3{ 1 };
	uint dim_t4{ 2 };
	uint cst_t4{ 1 };
	uint dim_t5{ 2 };
	uint cst_t5{ 1 };

	GKLSG_Int_Void __GetDim = (GKLSG_Int_Void)GetProcAddress(hDLL, "GetDim");
	GKLSG_Int_Void __GetCountCondition  = (GKLSG_Int_Void)GetProcAddress(hDLL, "GetCountCondition");

	uint dep{ 5 };
	Problem testProblem(2, 1, { 0.0, 0.0 }, { 1.0, 1.0 }, &f);
	Problem testProblem1(dim_t1, cst_t1, { 0.0, -1.0 }, { 4.0, 3.0 }, &task1);
	Problem testProblem2(dim_t2, cst_t2, { -2.5, -1.5 }, { 2.5, 1.5 }, &task2);
	Problem testProblem3(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task3);
	Problem testProblem4(dim_t4, cst_t4, { -1.0, -1.0 }, { 1.0, 1.0 }, &task4);
	Problem testProblem5(dim_t5, cst_t5, { -10, -10 }, { 10, 10 }, &task5);
	Problem testProblem6(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task6);

	//double share_solved_tasks = 0.0;
	//uint solved_tasks = 0;
	//double Delta = 1e-4;

	uint dim{ static_cast<uint>(__GetDim()) };
	uint cst{ static_cast<uint>(__GetCountCondition()) };

	//uint lim_iter{ 5000 };

	//std::vector<uint> iterations;
	//std::vector<int> flags(100);
	//for (auto& elem : flags) elem = 0;

	for (int i = 20; i < 21; ++i) {
		generate_task(i);
		std::cout << "\t\t\tTASK #" << i << std::endl;
		CoordinatesValues solution_point{ get_minimum_point() };
		double solution_minimum = global_minimum(solution_point);
		std::cout << "True global minimum is " << solution_minimum << std::endl;
		std::cout << "was reached at the point:" << std::endl;

		for (auto& elem : solution_point)
			std::cout << elem << std::endl;

		Problem testDll(dim, cst, { -1, -1 }, { 1, 1 }, &test_dll);

		uint max_iter{ 500 };
		double globalObj{ 4.0 };
		double globalCst{ 4.0 };
		double localObj{ 2.5 };
		double localCst{ 2.5 };
		double delta{ 1e-10 };
		double beta{ 0.4 };
		double eps{ 1e-9 };
		double diag{ static_cast<double>(MAX_POWER_THREE) };
		eps = eps * sqrt(static_cast<double>(dim)) * diag;

		Parameters param{ dim, cst, dep, max_iter, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1.0, 1.0 };
		std::shared_ptr<DivideByThree> solver(create_solver("SimplexMethodS", dim, cst, param, testDll));

		solver->solve();
		solver->write_generated_points();
		solver->write_generated_intervals();
		CoordinatesValues numerical_point{ solver->get_min_point() };

		std::cout << std::endl;
	}

	//std::ofstream out;
	//out.open("..\\x64\\Release\\NewTransformPM.txt");
	//if (out.is_open())
	//{
	//	for (auto& elem : iterations)
	//		out << elem << std::endl;
	//}

	//share_solved_tasks = 0.01 * iterations.size();
	//std::cout << share_solved_tasks << std::endl;

	/*uint dim{ dim_t3 };
	uint cst{ cst_t3 };
	uint max_iter{ 500 };
	double globalObj{ 4.5 };
	double globalCst{ 4.5 };
	double localObj{ 2.5 };
	double localCst{ 2.5 };
	double delta{ 1e-10 };
	double beta{ 0.4 };
	double eps{ 1e-9 };
	double diag{ static_cast<double>(MAX_POWER_THREE) };
	eps = eps * sqrt(static_cast<double>(dim)) * diag;

	Parameters param{ dim, cst, dep, max_iter, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1.5, 0.5, 1e-3, 100 };
	std::shared_ptr<DivideByThree> solver(create_solver("SimplexMethodS", dim, cst, param, testProblem3));

	solver->solve();
	solver->write_generated_points();
	solver->write_generated_intervals();*/

	deinit_dll();
	return EXIT_SUCCESS;
}