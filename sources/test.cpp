#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <Windows.h>

#include "synonymous_types.h"
#include "Functions.h"
#include "Problem.h"
#include "Parameters.h"
#include "DivideByThree.h"
#include "SolverFactory.h"

struct GKLSGParametersPacket {
	int numberProblem; // Было: targetFunction;
	int dimention; // Переставлено с первой позиции в списке !
	int maxNumberConstraints;    // Было: numberConstraints
	int numberLocalMinimasMinimandFunction;    // Было: numberLocalMinimasSourceFunction
	int numberLocalMinimasConstraintFunction;
	int minNumberActiveConstraints;
	int maxNumberActiveConstraints;
	double distanceFromGlobalMinimizerToParaboloidMinimizerMinimandFunction; //Новое поле !
	double minDistanceFromMinimizerToParaboloidMinimizer;
	double maxDistanceFromMinimizerToParaboloidMinimizer;
	double radiusAttractionRegionGlobalMinimizerMminimandFunction;  // Новое поле !
	double minRadiusAttractionRegionGlobalMinimizer;
	double maxRadiusAttractionRegionGlobalMinimizer;
	double valueGlobalMinimumMinimandFunction;    // Новое поле !
	double minValueGlobalMinimumConstraints;     // Было: minValueGlobalMinimum
	double maxValueGlobalMinimumConstraints;     // Было: maxValueGlobalMinimum
	double minAnglePhi;
	double maxAnglePhi;
	double maxPower;  // Было: maxP
	double probabilityConvexityConstraints;  // Было: beta
	double minNu0;   // Новое поле !
	double maxNu0;  // Было: nu0
	double minDeltaGlobal;  // Новое поле !
	double maxDeltaGlobal;  // Было: deltaGlobal
	double deltaLocal;  // Переставлено с последней позиции в списке !
	double maxLambda;
	double probabilityGglobalMinimumIsOutsidePermissibleRegion;  //Было:  globalBoundaryRate
	double relativeMeasurePermissibleRegion;    // Было:  measureAdmissibleRegion;
	int functionType;
};

struct SummaryPacket {
	wchar_t* summary;
	int length;
};

typedef int (WINAPIV* LPFN_Init) ();
typedef int (WINAPIV* LPFN_SetParam) (GKLSGParametersPacket*);
typedef int (WINAPIV* LPFN_Generate) (int);
typedef double (WINAPIV* LPFN_Solution) ();
typedef double (WINAPIV* LPFN_Obj) (double*);
typedef double (WINAPIV* LPFN_Cst) (double*, int);
typedef int (WINAPIV* LPFN_Free) ();
typedef int (WINAPIV* LPFN_Dim) ();
typedef int (WINAPIV* LPFN_Cond) ();
typedef int (WINAPIV* LPFN_Sum) (SummaryPacket* packet);

int main() {
	HMODULE hDLL = LoadLibrary(L"GKLSGLibrary.dll");
	if (hDLL == NULL) throw;

	GKLSGParametersPacket __parameters{
		1,
		2,
		4,
		10,
		20,
		2,
		2,
		0.9,
		0.2,
		1.0,
		0.2,
		0.1,
		0.3,
		-1.0,
		-2.9,
		-2.1,
		0.2,
		0.8,
		2.0,
		0.3,
		0.05,
		0.1,
		0.01,
		0.8,
		0.001,
		0.4,
		0.9,
		0.1,
		1
	};
	SummaryPacket __summary{ nullptr, 0 };

	LPFN_Init __Init				 = (LPFN_Init)GetProcAddress(hDLL, "Init");
	LPFN_SetParam   __SetGKLSGParameters = (LPFN_SetParam)GetProcAddress(hDLL, "SetGKLSGParameters");
	LPFN_Generate __GenerateGKLSG_Task = (LPFN_Generate)GetProcAddress(hDLL, "GenerateGKLSG_Task");
	LPFN_Free __FreeResources		 = (LPFN_Free)GetProcAddress(hDLL, "FreeResources");
	LPFN_Obj __TargetFunction	 = (LPFN_Obj)GetProcAddress(hDLL, "TargetFunction");
	LPFN_Cst __ConditionFunction  = (LPFN_Cst)GetProcAddress(hDLL, "ConditionFunction");
	LPFN_Dim __GetDim			 = (LPFN_Dim)GetProcAddress(hDLL, "GetDim");
	LPFN_Cond __GetCountCondition  = (LPFN_Cond)GetProcAddress(hDLL, "GetCountCondition");
	LPFN_Solution __GetGKLSGSolution = (LPFN_Solution)GetProcAddress(hDLL, "GetGKLSGSolution");
	LPFN_Sum __GetGKLSGSummary = (LPFN_Sum)GetProcAddress(hDLL, "GetGKLSGSummary");

	std::cout << __Init() << std::endl;
	__SetGKLSGParameters(&__parameters);
	int dm = __GetDim();
	int cd = __GetCountCondition();

	for (int i = 1; i < cd + 2; ++i)
		__GenerateGKLSG_Task(i);

	__GetGKLSGSummary(&__summary);
	std::cout << __GetGKLSGSolution() << std::endl;

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

	uint dep{ 5 };
	Problem testProblem(2, 1, { 0.0, 0.0 }, { 1.0, 1.0 }, &f);
	Problem testProblem1(dim_t1, cst_t1, { 0.0, -1.0 }, { 4.0, 3.0 }, &task1);
	Problem testProblem2(dim_t2, cst_t2, { -2.5, -1.5 }, { 2.5, 1.5 }, &task2);
	Problem testProblem3(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task3);
	Problem testProblem4(dim_t4, cst_t4, { -1.0, -1.0 }, { 1.0, 1.0 }, &task4);
	Problem testProblem5(dim_t5, cst_t5, { -10, -10 }, { 10, 10 }, &task5);
	Problem testProblem6(dim_t3, cst_t3, { 0, 0 }, { 2 * M_PI, 2 * M_PI }, &task6);

	uint dim{ dim_t1 };
	uint cst{ cst_t1 };
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
	std::shared_ptr<DivideByThree> solver(create_solver("LagrangeS", dim, cst, param, testProblem1));

	solver->solve();
	solver->write_generated_points();
	solver->write_generated_intervals();

	FreeLibrary(hDLL);

	return EXIT_SUCCESS;
}