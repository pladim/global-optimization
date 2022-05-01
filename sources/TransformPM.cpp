#include "TransformPM.h"

TransformPM::TransformPM(const uint& dimension,
						 const uint& constraints,
						 Parameters& parameters,
						 Problem& problem) :
	DivideByThree(dimension, constraints, parameters, problem),
	_areAllCharInfty(false),
	_doesGlobalChange(false),
	_localLipshEval(constraints + 1),
	_globalLipshEval(constraints + 1),
	_residual_minimum(std::numeric_limits<double>::max()) {}

void TransformPM::update_minimum(const FunctionsValues& evals,
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

double TransformPM::golden_ratio(double a, 
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

double TransformPM::calculate_residual(const double& t, const double& e, 
									   const uint& id_hyp) {
	Hyperinterval& hyp = _intervals[id_hyp];

	size_t idx_A = hyp.get_evalA();
	size_t idx_B = hyp.get_evalB();
	double mixed_LipshEval = 0.0;

	double result = -0.5 * (_evaluations[idx_A] - _current_minimum) * (t - e) / e;
	result += 0.5 * (_evaluations[idx_B] - _current_minimum) * (t + e) / e;
	mixed_LipshEval = mixedLipEval(hyp, 0);
	balance(mixed_LipshEval);
	result += mixed_LipshEval * (t * t - e * e);

	double candidate = 0.0;
	for (size_t i = 1; i < _constraints + 1; ++i) {
		candidate = -0.5 * _evaluations[idx_A + i] * (t - e) / e;
		candidate += 0.5 * _evaluations[idx_B + i] * (t + e) / e;
		mixed_LipshEval = mixedLipEval(hyp, i);
		balance(mixed_LipshEval);
		candidate += mixed_LipshEval * (t * t - e * e);

		if (result < candidate) result = candidate;
	}

	return result;
}

void TransformPM::calculate_characteristic(const uint& id_hyp) {
	if (_intervals[id_hyp].get_divisions() == 0) return;

	Hyperinterval& hyp = _intervals[id_hyp];
	uint index = (hyp.get_divisions() - 1) / _dimension;
	uint j = (hyp.get_divisions() - 1) % _dimension + 1;
	double h1 = static_cast<double>(HYPER_INTERVAL_SIDE_LENGTHS[20 - index]);
	double h2 = static_cast<double>(HYPER_INTERVAL_SIDE_LENGTHS[20 - index - 1]);
	double e = 0.5 * sqrt((_dimension - j) * h1 * h1 + j * h2 * h2);

	hyp.set_charact(golden_ratio(-e, e, e, id_hyp));
}