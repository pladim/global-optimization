#include <iostream>
#include "ParabolaMethod.h"

ParabolaMethod::ParabolaMethod(const uint& dimension,
							   const uint& constraints,
							   Parameters& parameters,
							   Problem& problem) :
	DivideByThree(dimension, constraints, parameters, problem),
	_areAllCharInfty(false),
	_doesGlobalChange(false),
	_localLipshEval(constraints + 1),
	_globalLipshEval(constraints + 1) {}

void ParabolaMethod::calculate_globalLipshConst(const uint& id_hyp) {
	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 2];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 1];

	for (uint i = 0; i < _constraints + 1; ++i) {
		if (_globalLipshEval[i] < hyp1.get_maxLipshEvaluations()[i]) {
			_doesGlobalChange = true;
			_globalLipshEval[i] = hyp1.get_maxLipshEvaluations()[i];
		}

		if (_globalLipshEval[i] < hyp2.get_maxLipshEvaluations()[i]) {
			_doesGlobalChange = true;
			_globalLipshEval[i] = hyp2.get_maxLipshEvaluations()[i];
		}

		if (_globalLipshEval[i] < hyp3.get_maxLipshEvaluations()[i]) {
			_doesGlobalChange = true;
			_globalLipshEval[i] = hyp3.get_maxLipshEvaluations()[i];
		}
	}

	if (_doesGlobalChange) {
		update_all_charact();
		_doesGlobalChange = false;
	}
	else {
		calculate_characteristic(id_hyp);
		calculate_characteristic(_generated_intervals - 2);
		calculate_characteristic(_generated_intervals - 1);
	}

	_areAllCharInfty = true;
	for (uint i = 0; i < _generated_intervals; ++i)
		if (_intervals[i].get_charact() < std::numeric_limits<double>::max())
			_areAllCharInfty = false;
}

double ParabolaMethod::mixedLipEval(const Hyperinterval& hyp, const uint& i) {
	double mixed_LipshitzEval = 0.0;
	double ratio = hyp.get_diagonal() / _parameters._criticalSize;
	double gained_global = 0.0;
	double gained_local = 0.0;
	double gain_global = 0.0;
	double gain_local = 0.0;

	if (i == 0) {
		gain_global = _parameters._gainGlobalObj;
		gain_local = _parameters._gainLocalObj;
	}
	else {
		gain_global = _parameters._gainGlobalCst;
		gain_local = _parameters._gainLocalCst;
	}

	if (_iteration >= _parameters._iter_thr) {
		gained_global = _globalLipshEval[i] + hyp.get_add_const(i);
		gained_local = std::min(
				gain_local * hyp.get_maxLipshEvaluations()[i],
				hyp.get_maxLipshEvaluations()[i] + hyp.get_add_const(i)
			);
	}
	else {
		gained_global = gain_global * _globalLipshEval[i];
		gained_local = gain_local * hyp.get_maxLipshEvaluations()[i];
	}

	if (hyp.get_diagonal() < _parameters._criticalSize) {
		mixed_LipshitzEval = ratio * gained_global + (1 - ratio) * gained_local;
	}
	else mixed_LipshitzEval = gained_global;

	return mixed_LipshitzEval;
}

void ParabolaMethod::update_all_charact() {
	for (uint i = 0; i < _generated_intervals; ++i)
		calculate_characteristic(i);
}

uint ParabolaMethod::min_by_charact() {
	auto it = std::min_element(
		_intervals.begin(),
		_intervals.end(),
		[&](const Hyperinterval& h1, const Hyperinterval& h2) {
			if ((h1.get_divisions() == MAX_EXPONENT_THREE * _dimension) || (h2.get_divisions() == MAX_EXPONENT_THREE * _dimension))
				return false;
			return h1.get_charact() < h2.get_charact();
		}
	);

	if (it == _intervals.end()) return _generated_intervals + 1;
	return it->get_id();
}

uint ParabolaMethod::max_by_length() {
	auto it = std::max_element(
		_intervals.begin(),
		_intervals.end(),
		[&](const Hyperinterval& h1, const Hyperinterval& h2) {
			if ((h1.get_divisions() == MAX_EXPONENT_THREE * _dimension) || (h2.get_divisions() == MAX_EXPONENT_THREE * _dimension))
				return false;
			return h1.get_diagonal() < h2.get_diagonal();
		}
	);

	if (it == _intervals.end()) return _generated_intervals + 1;
	return it->get_id();
}

uint ParabolaMethod::maxmin_lenchar() {
	auto it = std::max_element(
		_intervals.begin(),
		_intervals.end(),
		[](const Hyperinterval& h1, const Hyperinterval& h2) { return h1.get_diagonal() < h2.get_diagonal(); }
	);

	auto jt = std::min_element(
		_intervals.begin(),
		_intervals.end(),
		[&](const Hyperinterval& h1, const Hyperinterval& h2) { 
			if ((h1.get_divisions() == MAX_EXPONENT_THREE * _dimension) || (h2.get_divisions() == MAX_EXPONENT_THREE * _dimension))
				return false;

			if (std::abs(h1.get_diagonal() - it->get_diagonal()) < std::numeric_limits<double>::epsilon())
				return h1.get_charact() < it->get_charact();
			else return false;
		}
	);

	if (jt == _intervals.end()) return _generated_intervals + 1;
	return jt->get_id();
}

uint ParabolaMethod::supplementary() {
	std::vector<uint> cds;

	for (uint i = 0; i < _generated_intervals; ++i)
		if ((_intervals[i].get_idA() == _id_minimum) || (_intervals[i].get_idB() == _id_minimum))
			cds.push_back(i);

	auto it = std::min_element(
		cds.begin(),
		cds.end(),
		[&](const uint& i1, const uint& i2) {
			if ((_intervals[i1].get_divisions() == MAX_EXPONENT_THREE * _dimension) || (_intervals[i2].get_divisions() == MAX_EXPONENT_THREE * _dimension))
				return false;
			return _intervals[i1].get_charact() < _intervals[i2].get_charact(); 
		}
	);

	if (it == cds.end()) return _generated_intervals + 1;
	return *it;
}

uint ParabolaMethod::optimal_to_trisect() {
	if (_areAllCharInfty) return max_by_length();
	if (_iteration % 15 == 0) return maxmin_lenchar();
	if (_changes > 100) return supplementary();

	return min_by_charact();
}

uint ParabolaMethod::iterate(const uint& id_hyp) {
	if (_iteration == _parameters._iter_thr)
		for (uint i = 0; i < _generated_intervals; ++i)
			_intervals[i].count_add_const(_parameters._gainGlobalObj,
				                          _parameters._gainGlobalCst
			                             );
	trisect_interval(id_hyp);
	calculate_localLipshConst(id_hyp);
	calculate_globalLipshConst(id_hyp);
	return optimal_to_trisect();
}

void ParabolaMethod::balance(double& lipshConst) const {
	if (_iteration % 5 == 0) lipshConst *= _parameters._Gain;
	else lipshConst *= _parameters._Reduce;
}

void ParabolaMethod::solve() {
	initialization();
	uint id_current_interval = 0;
	bool flag = true;

	switch (_state) {
	case Mode::stop_by_precision:
		for (_iteration = 0; (_iteration < _max_it) && flag; ++_iteration) {
			id_current_interval = iterate(id_current_interval);
			if (id_current_interval == _generated_intervals + 1) break;
			Hyperinterval& hyp = _intervals[id_current_interval];

			if (_intervals[id_current_interval].get_diagonal() < _parameters._eps)
				flag = false;
		}

		if (!flag) {
			std::cout << "STOPPED BY PRECISION" << std::endl;
			_solved = true;
		}

		break;
	case Mode::test:
		bool adflag = true;

		for (_iteration = 0; ((_iteration < _max_it) && (adflag) && (flag)); ++_iteration) {
			id_current_interval = iterate(id_current_interval);
			if (id_current_interval == _generated_intervals + 1) break;
			Hyperinterval& hyp = _intervals[id_current_interval];

			if (_intervals[id_current_interval].get_diagonal() < _parameters._eps)
				adflag = false;

			if (is_close_enough()) {
				flag = false;
				_solved = true;
				std::cout << "SOLVED FOR TEST" << std::endl;
			}
		}

		if (!adflag) {
			std::cout << "STOPPED BY PRECISION" << std::endl;
		}

		break;
	}

	show_info();
}

bool ParabolaMethod::is_close_enough() {
	EncodedCoordinates ec(_dimension);
	Point& point = _points[_id_minimum];

	for (uint i = 0; i < _dimension; ++i)
		ec[i] = _coords[point.get_id_coord() + i];

	CoordinatesValues dc = _problem.decode_coordinates(ec);
	CoordinatesValues tm = _problem.get_true_minimum();

	double max_diff = std::fabs(-tm[0] + dc[0]);

	for (uint j = 1; j < _dimension; ++j) {
		double test_diff = std::fabs(-tm[j] + dc[j]);
		if (max_diff < test_diff) max_diff = test_diff;
	}

	if (max_diff < sqrt(_parameters._Delta) * 2) return true;
	else return false;
}

void ParabolaMethod::show_info() {
	std::cout << "Current minimum: " << _current_minimum << std::endl;

	EncodedCoordinates ec(_dimension);
	Point& point = _points[_id_minimum];

	for (uint i = 0; i < _dimension; ++i)
		ec[i] = _coords[point.get_id_coord() + i];

	CoordinatesValues dc = _problem.decode_coordinates(ec);

	std::cout << "Point: " << std::endl;
	for (uint i = 0; i < _dimension; ++i)
		std::cout << dc[i] << std::endl;

	std::cout << "Num of evaluations: " << _generated_points << std::endl;
}