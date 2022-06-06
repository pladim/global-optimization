#include <iostream>

#include "Tester.h"
#include "DivideByThree.h"

int main() {
	uint dep{ 3 };
	uint dim{ 4 };
	uint cst{ 4 };
	uint max_iter{ 500 };
	double globalObj{ 2.5 };
	double globalCst{ 2.5 };
	double localObj{ 1.2 };
	double localCst{ 1.2 };
	double delta{ 1e-10 };
	double beta{ 0.4 };
	double eps{ sqrt(static_cast<double>(dim)) + 0.1 };
	double diag{ static_cast<double>(MAX_POWER_THREE) };

	Parameters param{ dim, cst, dep, max_iter, localObj, localCst, globalObj, globalCst, delta, beta * diag, eps, 1.0, 1.0, 1e-3, 100, Mode::stop_by_precision };
	Tester test("SimplexMethodS", param);
	test.solve_task(105);

	return EXIT_SUCCESS;
}