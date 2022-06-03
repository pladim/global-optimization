#include "Simple.h"

Simple::Simple(const uint& dimension,
			   const uint& constraints,
			   Parameters& parameters,
			   Problem& problem) :
	Piyavskii(dimension, constraints, parameters, problem) {}

void Simple::calculate_characteristic(const uint& id_hyp) {
	Hyperinterval& hyp = _intervals[id_hyp];

	double mixed_LipshEval = mixedLipEval(hyp, 0);
	balance(mixed_LipshEval);

	double s_min = 0.5 * (_evaluations[hyp.get_evalA()] - _evaluations[hyp.get_evalB()]);
	s_min = s_min / (mixed_LipshEval * hyp.get_diagonal());
	s_min += 0.5;

	double s_m = 0.0;
	double s_p = 1.0;

	give_borders(s_m, s_p, hyp);

	if (s_m > s_p) {
		hyp.set_charact(std::numeric_limits<double>::max());
	}
	else if ((s_m < s_min) && (s_min < s_p)) {
		hyp.set_charact(
			0.5 * (_evaluations[hyp.get_evalA()] + _evaluations[hyp.get_evalB()]) -
			0.5 * mixed_LipshEval * hyp.get_diagonal()
		);
	}
	else if (s_min > s_p) {
		hyp.set_charact(
			std::min(
				_evaluations[hyp.get_evalA()] - mixed_LipshEval * hyp.get_diagonal() * s_m,
				_evaluations[hyp.get_evalA()] - mixed_LipshEval * hyp.get_diagonal() * s_p
			)
		);
	}
	else {
		hyp.set_charact(
			std::min(
				_evaluations[hyp.get_evalB()] - mixed_LipshEval * hyp.get_diagonal() * (1.0 - s_m),
				_evaluations[hyp.get_evalB()] - mixed_LipshEval * hyp.get_diagonal() * (1.0 - s_p)
			)
		);
	}
}

void Simple::give_borders(double& l, double& r, Hyperinterval& hyp) {
	double s_m = 0.0;
	double s_p = 1.0;
	double mixed_LipshEval = 0.0;

	for (uint i = 1; i < _constraints; ++i) {
		mixed_LipshEval = mixedLipEval(hyp, i);
		balance(mixed_LipshEval);

		s_m = _evaluations[hyp.get_evalA()] / (mixed_LipshEval * hyp.get_diagonal());
		s_p = 1 - (_evaluations[hyp.get_evalB()] / (mixed_LipshEval * hyp.get_diagonal()));

		if (s_m > s_p) {
			l = 1.0;
			r = 0.0;
			return;
		}
		else {
			if (s_m > l) l = s_m;
			if (s_p < r) r = s_p;
		}
	}
}

void Simple::update_minimum(const FunctionsValues& evals, 
							const uint& idp) {
	if (!_areAllCharInfty) {
		if (evals[0] < _current_minimum) {
			bool flag = true;
			for (uint i = 1; i < _constraints + 1; ++i)
				if (evals[i] > 0) flag = false;

			if (flag) {
				_current_minimum = evals[0];
				calc_distance(idp);
				_id_minimum = idp;
			}
		}
		else calc_distance(_id_minimum);
	}
	else {
		double candidate = *std::max_element(evals.begin() + 1, evals.end());

		if (candidate < _current_minimum) {
			_current_minimum = candidate;
			calc_distance(idp);
			_id_minimum = idp;
		}
		else calc_distance(_id_minimum);
	}
}