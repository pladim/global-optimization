#include <iostream>

#include "Tester.h"
#include "DivideByThree.h"

int main() {
	uint dep{ 3 };
	uint dim{ 1 };
	uint cst{ 1 };
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

	Parameters param{ dim, cst, dep, max_iter, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1.0, 1.0, 1e-3, 100, Mode::test };
	Tester test("SimplexMethodS", param);
	test.solve_task(105);
	test.write_measurements_to_file();

	return EXIT_SUCCESS;
}