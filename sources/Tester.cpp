#include "Tester.h"
#include "Functions.h"

Tester::Tester(const std::string& method, Parameters& param, const size_t& limit_iter, const double& Delta) :
	_limit_iter(limit_iter),
	_method(method),
	_param(param),
	_solution_point(param._dimension),
	_numerical_point(param._dimension),
	_measure_solved(limit_iter / 5),
	_Delta(Delta) {
	init_dll();
}

void Tester::run_test(const size_t& max_iter) {
	_param._max_it = max_iter;
	int solved_tasks = 0;
	for (int i = 1; i < 101; ++i) {
		generate_task(i);

		Problem pr(_param._dimension, _param._constraints, { -1, -1 }, { 1, 1 }, &test_dll);
		_problem = std::make_shared<Problem>(pr);
		_solver = create_solver(_method, _param._dimension, _param._constraints, _param, *_problem.get());

		_solver->solve();

		_solution_point = get_minimum_point();
		_numerical_point = _solver->get_min_point();

		double max_diff = std::fabs(-_solution_point[0] + _numerical_point[0]);
		for (uint j = 1; j < _param._dimension; ++j) {
			double test_diff = std::fabs(-_solution_point[j] + _numerical_point[j]);
			if (max_diff < test_diff) max_diff = test_diff;
		}

		if (max_diff < sqrt(_Delta) * 2) {
			++solved_tasks;
		}

		_solver.reset();
		_problem.reset();
	}

	_measure_solved[max_iter / 5 - 1] = 0.01 * solved_tasks;
}

void Tester::run() {
	for (size_t i = 5; i < _limit_iter; i += 5) {
		run_test(i);
	}
}

void Tester::write_results() {

}

Tester::~Tester() {
	deinit_dll();
}