#include <iostream>
#include "SimplePMwithSM.h"

SimplePMwithSM::SimplePMwithSM(const uint& dimension,
							   const uint& constraints,
							   Parameters& parameters,
							   Problem& problem) :
	WithSimplexMethod(dimension, constraints, parameters, problem) {}

void SimplePMwithSM::calculate_characteristic(const uint& id_hyp) {
	if (_intervals[id_hyp].get_divisions() == 0) return;
	Hyperinterval& hyp = _intervals[id_hyp];

	double mixed_LipshEval = mixedLipEval(hyp, 0);
	balance(mixed_LipshEval);
	double t_min = 0.5 * (_evaluations[hyp.get_evalA()] -
				   _evaluations[hyp.get_evalB()]);
	t_min = t_min / mixed_LipshEval;

	uint index = (hyp.get_divisions() - 1) / _dimension;
	uint j = (hyp.get_divisions() - 1) % _dimension + 1;
	double h1 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index];
	double h2 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index - 1];
	double e = 0.5 * sqrt((_dimension - j) * h1 * h1 + j * h2 * h2);

	t_min = t_min / e;
	double left = -e;
	double right = e;
	give_borders(left, right, hyp);

	if ((left == e) && (right == -e)) {
		hyp.set_charact(std::numeric_limits<double>::max());
	}
	else {
		if ((t_min > left) && (t_min < right)) {
			double charact = -_evaluations[hyp.get_evalA()] * (t_min - e);
			charact = charact + _evaluations[hyp.get_evalB()] * (t_min + e);
			charact = 0.5 * charact / e;
			charact = charact + 0.5 * mixed_LipshEval * (t_min * t_min - e * e);
			hyp.set_charact(charact);
		}
		else {
			double charact1 = -_evaluations[hyp.get_evalA()] * (left - e);
			charact1 = charact1 + _evaluations[hyp.get_evalB()] * (left + e);
			charact1 = 0.5 * charact1 / e;
			charact1 = charact1 + 0.5 * mixed_LipshEval * (left * left - e * e);

			double charact2 = -_evaluations[hyp.get_evalA()] * (right - e);
			charact2 = charact2 + _evaluations[hyp.get_evalB()] * (right + e);
			charact2 = 0.5 * charact2 / e;
			charact2 = charact2 + 0.5 * mixed_LipshEval * (right * right - e * e);

			hyp.set_charact(std::min(charact1, charact2));
		}
	}
}

void SimplePMwithSM::give_borders(double& l,
								  double& r,
								  Hyperinterval& hyp) {
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double fa = 0.0;
	double fb = 0.0;
	double M = 0.0;
	double e = r;

	for (uint i = 1; i < _constraints + 1; ++i) {
		M = mixedLipEval(hyp, i);
		fa = _evaluations[hyp.get_evalA() + i];
		fb = _evaluations[hyp.get_evalB() + i];

		a = 0.5 * M;
		b = 0.5 * (fb - fa) / e;
		c = 0.5 * (fb + fa - M * e * e);

		if (b * b - 4 * a * c >= 0) {
			fa = (-b - sqrt(b * b - 4 * a * c)) * 0.5 / a;
			fb = (-b + sqrt(b * b - 4 * a * c)) * 0.5 / a;

			if (fa > fb) std::swap(fa, fb);
			if ((fa >= l) && (fa <= r)) l = fa;
			if ((fb <= r) && (fb >= l)) r = fb;
		}
		else {
			l = e;
			r = -e;
			return;
		}
	}
}

void SimplePMwithSM::update_minimum(const FunctionsValues& evals,
									const uint& idp) {
	if (!_areAllCharInfty) {
		if (evals[0] < _current_minimum) {
			bool flag = true;
			for (uint i = 1; i < _constraints + 1; ++i)
				if (evals[i] > 0) flag = false;

			if (flag) {
				_current_minimum = evals[0];
				_id_minimum = idp;
			}
		}
	}
	else {
		double candidate = *std::max_element(evals.begin() + 1, evals.end());

		if (candidate < _current_minimum) {
			_current_minimum = candidate;
			_id_minimum = idp;
		}
	}
}