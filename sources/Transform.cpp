#include "Transform.h"

Transform::Transform(const uint& dimension,
					 const uint& constraints,
					 Parameters& parameters,
					 Problem& problem) :
	Piyavskii(dimension, constraints, parameters, problem) {}

void Transform::calculate_characteristic(const uint& id_hyp) {
	Hyperinterval& hyp = _intervals[id_hyp];
	hyp.set_charact(golden_ratio(0.0, 1.0, 1.0, id_hyp));
}

double Transform::golden_ratio(double a,
							   double b,
							   const double& e,
							   const uint& id_hyp) {
	if (b < a) std::swap(a, b);
	double x1 = b - phi * (b - a);
	double x2 = a + phi * (b - a);

	while (0.5 * (b - a) > 1e-10)
		if (calculate_residual(x1, e, id_hyp) > calculate_residual(x2, e, id_hyp)) {
			a = x1;
			x1 = b - (x2 - a);
		}
		else {
			b = x2;
			x2 = a + (b - x1);
		}

	return calculate_residual(0.5 * (b - a), e, id_hyp);
}

double Transform::calculate_residual(const double& t, const double& e,
									 const uint& id_hyp) {
	Hyperinterval& hyp = _intervals[id_hyp];
	double result = 0.0;

	size_t idx_A = hyp.get_evalA();
	size_t idx_B = hyp.get_evalB();
	double mixed_LipshEval = 0.0;
	mixed_LipshEval = mixedLipEval(hyp, 0);
	balance(mixed_LipshEval);

	double s_min = 0.5 * (_evaluations[hyp.get_evalA()] - _evaluations[hyp.get_evalB()]);
	s_min = s_min / (mixed_LipshEval * hyp.get_diagonal());
	s_min += 0.5;

	if (t < s_min) {
		result = _evaluations[hyp.get_evalA()] - _current_minimum - mixed_LipshEval * hyp.get_diagonal() * t;
	}
	else {
		result = _evaluations[hyp.get_evalB()] - _current_minimum - mixed_LipshEval * hyp.get_diagonal() * (1 - t);
	}

	double candidate = 0.0;
	for (size_t i = 1; i < _constraints + 1; ++i) {
		mixed_LipshEval = mixedLipEval(hyp, i);
		balance(mixed_LipshEval);

		s_min = 0.5 * (_evaluations[hyp.get_evalA() + i] - _evaluations[hyp.get_evalB() + i]);
		s_min = s_min / (mixed_LipshEval * hyp.get_diagonal());
		s_min += 0.5;

		if (t < s_min) {
			candidate = _evaluations[hyp.get_evalA() + i] - mixed_LipshEval * hyp.get_diagonal() * t;
		}
		else {
			candidate = _evaluations[hyp.get_evalB() + i] - mixed_LipshEval * hyp.get_diagonal() * (1 - t);
		}

		if (result < candidate) result = candidate;
	}

	return result;
}

void Transform::update_minimum(const FunctionsValues& evals,
	const uint& idp) {
	if (evals[0] < _current_minimum) {
		bool flag = true;
		for (uint i = 1; i < _constraints + 1; ++i)
			if (evals[i] > 0) flag = false;

		if (flag) {
			_current_minimum = evals[0];
			_id_minimum = idp;
			update_all_charact();
		}
	}
}