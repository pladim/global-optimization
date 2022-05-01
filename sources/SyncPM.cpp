//#include <iostream>
//#include "SyncPM.h"
//
//SyncPM::SyncPM(const uint& dimension,
//			   const uint& constraints,
//			   Parameters& parameters,
//			   Problem& problem) :
//	DivideByThree(dimension, constraints, parameters, problem),
//	_areAllCharInfty(false),
//	_doesGlobalChange(false),
//	_localLipshEval(constraints + 1),
//	_globalLipshEval(constraints + 1),
//	_point(4 * _dimension),
//	_incs(6 * 2),
//	_fval(4),
//	_non_proj_incs(6 * _dimension) {}
//
//double SyncPM::calculate_residual(const double& t, const double& e,
//							      const uint& id_hyp) {
//	Hyperinterval& hyp = _intervals[id_hyp];
//
//	size_t idx_A = hyp.get_evalA();
//	size_t idx_B = hyp.get_evalB();
//	double mixed_LipshEval = 0.0;
//
//	double result = -0.5 * (_evaluations[idx_A] - _current_minimum) * (t - e) / e;
//	result += 0.5 * (_evaluations[idx_B] - _current_minimum) * (t + e) / e;
//	mixed_LipshEval = mixedLipEval(hyp, 0);
//	balance(mixed_LipshEval);
//	result += mixed_LipshEval * (t * t - e * e);
//
//	double candidate = 0.0;
//	for (size_t i = 1; i < _constraints + 1; ++i) {
//		candidate = -0.5 * _evaluations[idx_A + i] * (t - e) / e;
//		candidate += 0.5 * _evaluations[idx_B + i] * (t + e) / e;
//		mixed_LipshEval = mixedLipEval(hyp, i);
//		balance(mixed_LipshEval);
//		candidate += mixed_LipshEval * (t * t - e * e);
//
//		if (result < candidate) result = candidate;
//	}
//
//	return result;
//}
//
//double SyncPM::golden_ratio(double a,
//							double b,
//							const double& e,
//							const uint& id_hyp) {
//	if (b < a) std::swap(a, b);
//	double x1 = b - phi * (b - a);
//	double x2 = a + phi * (b - a);
//
//	while (0.5 * (b - a) > 1e-10)
//		if (calculate_residual(x1, e, id_hyp) > calculate_residual(x2, e, id_hyp)) {
//			a = x1;
//			x1 = b - (x2 - a);
//		}
//		else {
//			b = x2;
//			x2 = a + (b - x1);
//		}
//
//	return calculate_residual(0.5 * (b - a), e, id_hyp);
//}
//
//void SyncPM::calculate_characteristic(const uint & id_hyp) {
//	if (_intervals[id_hyp].get_divisions() == 0) return;
//	Hyperinterval& hyp = _intervals[id_hyp];
//
//	double mixed_LipshEval = mixedLipEval(hyp, 0);
//	balance(mixed_LipshEval);
//	double t_min = 0.5 * (_evaluations[hyp.get_evalA()] -
//		_evaluations[hyp.get_evalB()]);
//	t_min = t_min / mixed_LipshEval;
//
//	uint index = (hyp.get_divisions() - 1) / _dimension;
//	uint j = (hyp.get_divisions() - 1) % _dimension + 1;
//	double h1 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index];
//	double h2 = (double)HYPER_INTERVAL_SIDE_LENGTHS[20 - index - 1];
//	double e = 0.5 * sqrt((_dimension - j) * h1 * h1 + j * h2 * h2);
//
//	t_min = t_min / e;
//	double left = -e;
//	double right = e;
//	give_borders(left, right, hyp);
//
//	if ((left == e) && (right == -e)) {
//		hyp.set_charact(std::numeric_limits<double>::max());
//	}
//	else {
//		if ((t_min > left) && (t_min < right)) {
//			double charact = -_evaluations[hyp.get_evalA()] * (t_min - e);
//			charact = charact + _evaluations[hyp.get_evalB()] * (t_min + e);
//			charact = 0.5 * charact / e;
//			charact = charact + 0.5 * mixed_LipshEval * (t_min * t_min - e * e);
//
//			double charact_ = golden_ratio(-e, e, e, id_hyp);
//			hyp.set_charact(std::min(charact, charact_));
//		}
//		else {
//			double charact1 = -_evaluations[hyp.get_evalA()] * (left - e);
//			charact1 = charact1 + _evaluations[hyp.get_evalB()] * (left + e);
//			charact1 = 0.5 * charact1 / e;
//			charact1 = charact1 + 0.5 * mixed_LipshEval * (left * left - e * e);
//
//			double charact2 = -_evaluations[hyp.get_evalA()] * (right - e);
//			charact2 = charact2 + _evaluations[hyp.get_evalB()] * (right + e);
//			charact2 = 0.5 * charact2 / e;
//			charact2 = charact2 + 0.5 * mixed_LipshEval * (right * right - e * e);
//
//			charact1 = std::min(charact1, charact2);
//			double charact_ = golden_ratio(-e, e, e, id_hyp);
//			hyp.set_charact(std::min(charact1, charact_));
//		}
//	}
//}
//
//void SyncPM::calculate_localLipshConst(const uint& id_hyp) {
//	Hyperinterval& hyp1 = _intervals[id_hyp];
//	Hyperinterval& hyp2 = _intervals[_generated_intervals - 2];
//	Hyperinterval& hyp3 = _intervals[_generated_intervals - 1];
//
//	uint pos_a = hyp1.get_idA();
//	uint pos_v = hyp2.get_idA();
//	uint pos_u = hyp2.get_idB();
//	uint pos_b = hyp3.get_idB();
//
//	decode_and_save(pos_a, 0);
//	decode_and_save(pos_v, 1);
//	decode_and_save(pos_u, 2);
//	decode_and_save(pos_b, 3);
//
//	calculate_and_project(hyp1.get_previous_axis());
//
//	for (uint i = 0; i < _constraints + 1; ++i) {
//		generate_right_part(i, pos_a * (_constraints + 1),
//			pos_v * (_constraints + 1),
//			pos_u * (_constraints + 1),
//			pos_b * (_constraints + 1));
//
//		get_solution();
//
//		_localLipshEval[i] = get_solution() / (static_cast<double>(MAX_POWER_THREE) * static_cast<double>(MAX_POWER_THREE));
//	}
//
//	hyp1.update_localLipQueues(_localLipshEval,
//		_parameters._delta /
//		((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
//	hyp2.update_localLipQueues(_localLipshEval,
//		_parameters._delta /
//		((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
//	hyp3.update_localLipQueues(_localLipshEval,
//		_parameters._delta /
//		((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
//}
//
//void SyncPM::decode_and_save(const uint& pos, const uint& order) {
//	for (size_t i = 0; i < _dimension; ++i)
//		_transit1[i] = _coords[pos * _dimension + i];
//	CoordinatesValues& t = _problem.decode_coordinates01(_transit1);
//	for (size_t i = 0; i < _dimension; ++i)
//		_point[i + _dimension * order] = t[i];
//}
//
//void SyncPM::projection(const uint& order,
//								std::vector<double>& e1,
//								std::vector<double>& e2) {
//	_incs[2 * order] = scalar_product(_non_proj_incs.begin() + order * _dimension,
//		e2.begin());
//	_incs[2 * order + 1] = scalar_product(_non_proj_incs.begin() + order * _dimension,
//		e1.begin());
//}
//
//void SyncPM::calculate_and_project(const uint& axis) {
//	// высчитываем приращения
//	for (uint i = 0; i < _dimension; ++i) {
//		// 21
//		_non_proj_incs[i] = _point[i + _dimension] - _point[i];
//		// 31
//		_non_proj_incs[i + _dimension] =
//			_point[i + 2 * _dimension] - _point[i];
//		// 41
//		_non_proj_incs[i + 2 * _dimension] =
//			_point[i + 3 * _dimension] - _point[i];
//		// 32
//		_non_proj_incs[i + 3 * _dimension] =
//			_point[i + 2 * _dimension] - _point[i + _dimension];
//		// 42
//		_non_proj_incs[i + 4 * _dimension] =
//			_point[i + 3 * _dimension] - _point[i + _dimension];
//		// 43
//		_non_proj_incs[i + 5 * _dimension] =
//			_point[i + 3 * _dimension] - _point[i + 2 * _dimension];
//	}
//
//	std::vector<double> e2(_dimension, 0);
//	e2[axis] = 1;
//	std::vector<double> e1(_dimension, 0);
//	double scalar = _non_proj_incs[2 * _dimension + axis];
//	for (uint i = 0; i < _dimension; ++i)
//		e1[i] = _non_proj_incs[i + 2 * _dimension] - scalar * e2[i];
//
//	scalar = scalar_product(e1.begin(), e1.begin());
//	for (uint i = 0; i < _dimension; ++i)
//		e1[i] /= sqrt(scalar);
//
//	projection(0, e1, e2);
//	projection(1, e1, e2);
//	projection(2, e1, e2);
//	projection(3, e1, e2);
//	projection(4, e1, e2);
//	projection(5, e1, e2);
//}
//
//double SyncPM::scalar_product(const std::vector<double>::iterator& a,
//	const std::vector<double>::iterator& b) {
//	double result = 0.0;
//
//	for (size_t i = 0; i < _dimension; ++i)
//		result += *(a + i) * *(b + i);
//
//	return result;
//}
//
//void SyncPM::generate_right_part(const uint& function,
//	const uint& eval_a,
//	const uint& eval_v,
//	const uint& eval_u,
//	const uint& eval_b) {
//	_fval[0] = _evaluations[eval_a + function];
//	_fval[1] = _evaluations[eval_v + function];
//	_fval[2] = _evaluations[eval_u + function];
//	_fval[3] = _evaluations[eval_b + function];
//}
//
//double SyncPM::get_solution() {
//	double result = 0.0;
//	Py_Initialize();
//	PyObject* _module = PyImport_ImportModule("solver");
//	PyObject* _function = PyObject_GetAttrString(_module, (char*)"solve");
//	PyObject* _args = PyTuple_Pack(16,
//		PyFloat_FromDouble(_incs[0]),
//		PyFloat_FromDouble(_incs[1]),
//		PyFloat_FromDouble(_incs[2]),
//		PyFloat_FromDouble(_incs[3]),
//		PyFloat_FromDouble(_incs[4]),
//		PyFloat_FromDouble(_incs[5]),
//		PyFloat_FromDouble(_incs[6]),
//		PyFloat_FromDouble(_incs[7]),
//		PyFloat_FromDouble(_incs[8]),
//		PyFloat_FromDouble(_incs[9]),
//		PyFloat_FromDouble(_incs[10]),
//		PyFloat_FromDouble(_incs[11]),
//		PyFloat_FromDouble(_fval[0]),
//		PyFloat_FromDouble(_fval[1]),
//		PyFloat_FromDouble(_fval[2]),
//		PyFloat_FromDouble(_fval[3])
//	);
//
//	PyObject* _result = PyObject_CallObject(_function, _args);
//	result = PyFloat_AsDouble(_result);
//
//	Py_Finalize();
//	return result;
//}
//
//double SyncPM::mixedLipEval(const Hyperinterval& hyp, const uint& i) {
//	double mixed_LipshitzEval = 0.0;
//	double ratio = hyp.get_diagonal() / _parameters._criticalSize;
//
//	if (i == 0) {
//		if (hyp.get_diagonal() < _parameters._criticalSize) {
//			mixed_LipshitzEval = ratio * _parameters._gainGlobalObj *
//				_globalLipshEval[i];
//			mixed_LipshitzEval += (1 - ratio) * _parameters._gainLocalObj *
//				hyp.get_maxLipshEvaluations()[i];
//		}
//		else mixed_LipshitzEval = _parameters._gainGlobalObj *
//			_globalLipshEval[i];
//	}
//	else {
//		if (hyp.get_diagonal() < _parameters._criticalSize) {
//			mixed_LipshitzEval = ratio * _parameters._gainGlobalCst *
//				_globalLipshEval[i];
//			mixed_LipshitzEval += (1 - ratio) * _parameters._gainLocalCst *
//				hyp.get_maxLipshEvaluations()[i];
//		}
//		else mixed_LipshitzEval = _parameters._gainGlobalCst *
//			_globalLipshEval[i];
//	}
//
//	return mixed_LipshitzEval;
//}
//
//void SyncPM::give_borders(double& l,
//						  double& r,
//						  Hyperinterval& hyp) {
//	double a = 0.0;
//	double b = 0.0;
//	double c = 0.0;
//	double fa = 0.0;
//	double fb = 0.0;
//	double M = 0.0;
//	double e = r;
//
//	for (uint i = 1; i < _constraints + 1; ++i) {
//		M = mixedLipEval(hyp, i);
//		fa = _evaluations[hyp.get_evalA() + i];
//		fb = _evaluations[hyp.get_evalB() + i];
//
//		a = 0.5 * M;
//		b = 0.5 * (fb - fa) / e;
//		c = 0.5 * (fb + fa - M * e * e);
//
//		if (b * b - 4 * a * c >= 0) {
//			fa = (-b - sqrt(b * b - 4 * a * c)) * 0.5 / a;
//			fb = (-b + sqrt(b * b - 4 * a * c)) * 0.5 / a;
//
//			if (fa > fb) std::swap(fa, fb);
//			if ((fa >= l) && (fa <= r)) l = fa;
//			if ((fb <= r) && (fb >= l)) r = fb;
//		}
//		else {
//			l = e;
//			r = -e;
//			return;
//		}
//	}
//}
//
//void SyncPM::balance(double& lipshConst) const {
//	if (_iteration % 5 == 0) lipshConst *= _parameters._Gain;
//	else lipshConst *= _parameters._Reduce;
//}
//
//void SyncPM::calculate_globalLipshConst(const uint& id_hyp) {
//	Hyperinterval& hyp1 = _intervals[id_hyp];
//	Hyperinterval& hyp2 = _intervals[_generated_intervals - 2];
//	Hyperinterval& hyp3 = _intervals[_generated_intervals - 1];
//
//	for (uint i = 0; i < _constraints + 1; ++i) {
//		if (_globalLipshEval[i] < hyp1.get_maxLipshEvaluations()[i]) {
//			_doesGlobalChange = true;
//			_globalLipshEval[i] = hyp1.get_maxLipshEvaluations()[i];
//		}
//
//		if (_globalLipshEval[i] < hyp2.get_maxLipshEvaluations()[i]) {
//			_doesGlobalChange = true;
//			_globalLipshEval[i] = hyp2.get_maxLipshEvaluations()[i];
//		}
//
//		if (_globalLipshEval[i] < hyp3.get_maxLipshEvaluations()[i]) {
//			_doesGlobalChange = true;
//			_globalLipshEval[i] = hyp3.get_maxLipshEvaluations()[i];
//		}
//	}
//
//	if (_doesGlobalChange) {
//		update_all_charact();
//		_doesGlobalChange = false;
//	}
//	else {
//		calculate_characteristic(id_hyp);
//		calculate_characteristic(_generated_intervals - 2);
//		calculate_characteristic(_generated_intervals - 1);
//	}
//}
//
//void SyncPM::update_all_charact() {
//	for (uint i = 0; i < _generated_intervals; ++i)
//		calculate_characteristic(i);
//}
//
//uint SyncPM::optimal_to_trisect() {
//	uint id_optimal_hyp = 0;
//
//	if (_areAllCharInfty) {
//		double optimal_charact = _intervals[id_optimal_hyp].get_diagonal();
//		double current_charact = 0.0;
//
//		for (uint id_hyp = 1; id_hyp < _generated_intervals; ++id_hyp) {
//			current_charact = _intervals[id_hyp].get_diagonal();
//			if (std::abs(current_charact - optimal_charact) < std::numeric_limits<double>::epsilon());
//			else if (optimal_charact < current_charact) {
//				optimal_charact = current_charact;
//				id_optimal_hyp = id_hyp;
//			}
//		}
//	}
//	else {
//		double optimal_charact = _intervals[id_optimal_hyp].get_charact();
//		double current_charact = 0.0;
//
//		for (uint id_hyp = 1; id_hyp < _generated_intervals; ++id_hyp) {
//			current_charact = _intervals[id_hyp].get_charact();
//			if (std::abs(current_charact - optimal_charact) < std::numeric_limits<double>::epsilon());
//			else if (optimal_charact > current_charact) {
//				optimal_charact = current_charact;
//				id_optimal_hyp = id_hyp;
//			}
//		}
//	}
//
//	return id_optimal_hyp;
//}
//
//uint SyncPM::iterate(const uint& id_hyp) {
//	trisect_interval(id_hyp);
//	calculate_localLipshConst(id_hyp);
//	calculate_globalLipshConst(id_hyp);
//	return optimal_to_trisect();
//}
//
//void SyncPM::solve() {
//	initialization();
//	uint id_current_interval = 0;
//	bool flag = true;
//
//	for (uint i = 0; ((i < 250) && (flag)); ++i) {
//		id_current_interval = iterate(id_current_interval);
//		Hyperinterval& hyp = _intervals[id_current_interval];
//
//		if (_intervals[id_current_interval].get_diagonal() < _parameters._eps)
//			flag = false;
//	}
//
//	if (!flag) std::cout << "STOPPED BY PRECISION" << std::endl;
//
//	std::cout << "Current minimum: " << _current_minimum << std::endl;
//	EncodedCoordinates ec(_dimension);
//	Point& point = _points[_id_minimum];
//	for (uint i = 0; i < _dimension; ++i)
//		ec[i] = _coords[point.get_id_coord() + i];
//	CoordinatesValues dc = _problem.decode_coordinates(ec);
//	std::cout << "Point: " << std::endl;
//	for (uint i = 0; i < _dimension; ++i)
//		std::cout << dc[i] << std::endl;
//	std::cout << "Num of evaluations: " << _generated_points << std::endl;
//}