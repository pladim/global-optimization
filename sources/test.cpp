#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

#include "synonymous_types.h"
#include "Functions.h"
#include "Problem.h"
#include "Parameters.h"
#include "DivideByThree.h"
#include "SolverFactory.h"

int main() {
	//LPFN_Generate __GenerateGKLSG_Task = (LPFN_Generate)GetProcAddress(hDLL, "GenerateGKLSG_Task");
	//LPFN_Free __FreeResources		 = (LPFN_Free)GetProcAddress(hDLL, "FreeResources");
	// LPFN_Obj __TargetFunction	 = (LPFN_Obj)GetProcAddress(hDLL, "TargetFunction");
	// LPFN_Cst __ConditionFunction  = (LPFN_Cst)GetProcAddress(hDLL, "ConditionFunction");
	// LPFN_Dim __GetDim			 = (LPFN_Dim)GetProcAddress(hDLL, "GetDim");
	// LPFN_Cond __GetCountCondition  = (LPFN_Cond)GetProcAddress(hDLL, "GetCountCondition");
	// LPFN_Solution __GetGKLSGSolution = (LPFN_Solution)GetProcAddress(hDLL, "GetGKLSGSolution");

	init_dll();

	// __GetGKLSGSummary(&__summary);

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

	LPFN_Dim __GetDim = (LPFN_Dim)GetProcAddress(hDLL, "GetDim");
	LPFN_Cond __GetCountCondition  = (LPFN_Cond)GetProcAddress(hDLL, "GetCountCondition");

	uint dep{ 5 };
	Problem testProblem(2, 1, { 0.0, 0.0 }, { 1.0, 1.0 }, &f);
	Problem testProblem1(dim_t1, cst_t1, { 0.0, -1.0 }, { 4.0, 3.0 }, &task1);
	Problem testProblem2(dim_t2, cst_t2, { -2.5, -1.5 }, { 2.5, 1.5 }, &task2);
	Problem testProblem3(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task3);
	Problem testProblem4(dim_t4, cst_t4, { -1.0, -1.0 }, { 1.0, 1.0 }, &task4);
	Problem testProblem5(dim_t5, cst_t5, { -10, -10 }, { 10, 10 }, &task5);
	Problem testProblem6(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task6);
	Problem testDll(__GetDim(), __GetCountCondition(), {-1, -1}, {1, 1}, &test_dll);

	uint dim{ static_cast<uint>(__GetDim())};
	uint cst{ static_cast<uint>(__GetCountCondition())};
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

	/*std::ofstream _log;
	_log.open("D:\\materials\\projects\\visualize_hyperinterval\\status.txt", std::ios::app);
	if (_log.is_open()) {
		while (globalObj < 5.0) {
			while (globalCst < 5.0) {
				while (localObj < 3.0) {
					while (localCst < 3.0) {
						_log << "TASK\t" << globalObj << " " << globalCst << " " << localObj << " " << localCst << endl;
						Parameters param{ dim, cst, dep, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1.0, 1.0 };
						std::shared_ptr<DivideByThree> solver(create_solver("SimplexMethodT", dim, cst, param, testProblem1));
						solver->solve();
						_log << solver->get_min() << std::endl;
						cout << "\n\n";

						localCst += 0.5;
					}

					localCst = 1.5;
					localObj += 0.5;
				}

				localCst = 1.5;
				localObj = 1.5;
				globalCst += 0.5;
			}

			localCst = 1.5;
			localObj = 1.5;
			globalCst = 2.5;
			globalObj += 0.5;
		}
	}*/

	Parameters param{ dim, cst, dep, max_iter, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1.0, 1.0 };
	std::shared_ptr<DivideByThree> solver(create_solver("SimplexMethodS", dim, cst, param, testDll));

	solver->solve();
	solver->write_generated_points();
	solver->write_generated_intervals();

	deinit_dll();
	return EXIT_SUCCESS;
}