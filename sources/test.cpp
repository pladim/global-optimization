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
#include "Tester.h"

int main() {
	init_dll();

	GKLSG_Int_Void __GetDim = (GKLSG_Int_Void)GetProcAddress(hDLL, "GetDim");
	GKLSG_Int_Void __GetCountCondition  = (GKLSG_Int_Void)GetProcAddress(hDLL, "GetCountCondition");

	uint dep{ 3 };
	Problem testProblem(2, 1, { 0.0, 0.0 }, { 1.0, 1.0 }, &f);

	uint dim{ static_cast<uint>(__GetDim()) };
	uint cst{ static_cast<uint>(__GetCountCondition()) };

	std::vector<uint> iterations;

	for (int i = 1; i < 101; ++i) {
		generate_task(i);
		std::cout << "\t\t\tTASK #" << i << std::endl;
		CoordinatesValues solution_point{ get_minimum_point() };
		double solution_minimum = global_minimum(solution_point);
		std::cout << "True global minimum is " << solution_minimum << std::endl;
		std::cout << "was reached at the point:" << std::endl;

		for (auto& elem : solution_point)
			std::cout << elem << std::endl;

		Problem testDll(dim, cst, { -1, -1 }, { 1, 1 }, &test_dll);
		testDll.set_true_minimum(solution_point);

		uint max_iter{ 1500 };
		double globalObj{ 3.0 };
		double globalCst{ 3.0 };
		double localObj{ 1.5 };
		double localCst{ 1.5 };
		double delta{ 1e-10 };
		double beta{ 0.4 };
		double eps{ 1e-9 };
		double diag{ static_cast<double>(MAX_POWER_THREE) };
		eps = eps * sqrt(static_cast<double>(dim)) * diag;

		Parameters param{ dim, cst, dep, max_iter, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1.0, 1.0, 1e-3, 100 };
		std::shared_ptr<DivideByThree> solver(create_solver("PiyavskiiT", dim, cst, param, testDll));

		solver->solve();
		if (solver->solved()) iterations.push_back(solver->get_gen());

		std::cout << std::endl;
	}

	/*uint dim{ dim_t3 };
	uint cst{ cst_t3 };
	uint max_iter{ 500 };
	double globalObj{ 2.5 };
	double globalCst{ 2.5 };
	double localObj{ 1.5 };
	double localCst{ 1.5 };
	double delta{ 1e-10 };
	double beta{ 0.4 };
	double eps{ 1e-9 };
	double diag{ static_cast<double>(MAX_POWER_THREE) };
	eps = eps * sqrt(static_cast<double>(dim)) * diag;

	Parameters param{ dim, cst, dep, max_iter, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1.0, 1.0, 1e-3, 100 };
	std::shared_ptr<DivideByThree> solver(create_solver("SimplexMethodT", dim, cst, param, testProblem3));

	solver->solve();
	solver->write_generated_points();
	solver->write_generated_intervals();*/

	deinit_dll();
	return EXIT_SUCCESS;
}