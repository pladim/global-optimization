#include <iostream>
#include "SimplePMwithSM.h"

SimplePMwithSM::SimplePMwithSM(const uint& dimension,
							   const uint& constraints,
							   Parameters& parameters,
							   Problem& problem) :
	SimplePM(dimension, constraints, parameters, problem),
	_point(4 * _dimension),
	_incs(6 * 2),
	_fval(4),
	_non_proj_incs(6 * _dimension) {
}

void SimplePMwithSM::decode_and_save(const uint& pos, const uint& order) {
	for (size_t i = 0; i < _dimension; ++i)
		_transit1[i] = _coords[pos * _dimension + i];
	CoordinatesValues& t = _problem.decode_coordinates01(_transit1);
	for (size_t i = 0; i < _dimension; ++i)
		_point[i + _dimension * order] = t[i];
}

void SimplePMwithSM::projection(const uint& order, 
								std::vector<double>& e1, 
								std::vector<double>& e2) {
	_incs[2 * order]     = scalar_product(_non_proj_incs.begin() + order * _dimension, 
									      e2.begin());
	_incs[2 * order + 1] = scalar_product(_non_proj_incs.begin() + order * _dimension, 
										  e1.begin());
}

void SimplePMwithSM::calculate_localLipshConst(const uint& id_hyp) {
	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 2];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 1];

	uint pos_a = hyp1.get_idA();
	uint pos_v = hyp2.get_idA();
	uint pos_u = hyp2.get_idB();
	uint pos_b = hyp3.get_idB();

	decode_and_save(pos_a, 0);
	decode_and_save(pos_v, 1);
	decode_and_save(pos_u, 2);
	decode_and_save(pos_b, 3);

	calculate_and_project(hyp1.get_previous_axis());
	
	for (uint i = 0; i < _constraints + 1; ++i) {
		generate_right_part(i, pos_a * (_constraints + 1),
							   pos_v * (_constraints + 1),
							   pos_u * (_constraints + 1),
							   pos_b * (_constraints + 1));

		get_solution();

		_localLipshEval[i] = get_solution() / (static_cast<double>(MAX_POWER_THREE) * static_cast<double>(MAX_POWER_THREE));
	}

	hyp1.update_localLipQueues(_localLipshEval, 
							   _parameters._delta / 
							   ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
	hyp2.update_localLipQueues(_localLipshEval, 
							   _parameters._delta / 
							   ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
	hyp3.update_localLipQueues(_localLipshEval, 
							   _parameters._delta / 
							   ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE));
}

void SimplePMwithSM::calculate_and_project(const uint& axis) {
	// высчитываем приращения
	for (uint i = 0; i < _dimension; ++i) {
		// 21
		_non_proj_incs[i] = _point[i + _dimension] - _point[i];
		// 31
		_non_proj_incs[i + _dimension] =
			_point[i + 2 * _dimension] - _point[i];
		// 41
		_non_proj_incs[i + 2 * _dimension] =
			_point[i + 3 * _dimension] - _point[i];
		// 32
		_non_proj_incs[i + 3 * _dimension] =
			_point[i + 2 * _dimension] - _point[i + _dimension];
		// 42
		_non_proj_incs[i + 4 * _dimension] =
			_point[i + 3 * _dimension] - _point[i + _dimension];
		// 43
		_non_proj_incs[i + 5 * _dimension] =
			_point[i + 3 * _dimension] - _point[i + 2 * _dimension];
	}

	std::vector<double> e2(_dimension, 0);
	e2[axis] = 1;
	std::vector<double> e1(_dimension, 0);
	double scalar = _non_proj_incs[2 * _dimension + axis];
	for (uint i = 0; i < _dimension; ++i)
		e1[i] = _non_proj_incs[i + 2 * _dimension] - scalar * e2[i];

	scalar = scalar_product(e1.begin(), e1.begin());
	for (uint i = 0; i < _dimension; ++i)
		e1[i] /= sqrt(scalar);

	projection(0, e1, e2);
	projection(1, e1, e2);
	projection(2, e1, e2);
	projection(3, e1, e2);
	projection(4, e1, e2);
	projection(5, e1, e2);
}

double SimplePMwithSM::scalar_product(const std::vector<double>::iterator& a,
									  const std::vector<double>::iterator& b) {
	double result = 0.0;

	for (size_t i = 0; i < _dimension; ++i)
		result += *(a + i) * *(b + i);

	return result;
}

void SimplePMwithSM::generate_right_part(const uint& function,
										 const uint& eval_a,
										 const uint& eval_v,
										 const uint& eval_u,
										 const uint& eval_b) {
	_fval[0] = _evaluations[eval_a + function];
	_fval[1] = _evaluations[eval_v + function];
	_fval[2] = _evaluations[eval_u + function];
	_fval[3] = _evaluations[eval_b + function];
}

double SimplePMwithSM::get_solution() {
	double result = 0.0;
	Py_Initialize();
	PyObject* _module = PyImport_ImportModule("solver");
	PyObject* _function = PyObject_GetAttrString(_module, (char*)"solve");
	PyObject* _args = PyTuple_Pack(16,
		PyFloat_FromDouble(_incs[0]),
		PyFloat_FromDouble(_incs[1]),
		PyFloat_FromDouble(_incs[2]),
		PyFloat_FromDouble(_incs[3]),
		PyFloat_FromDouble(_incs[4]),
		PyFloat_FromDouble(_incs[5]),
		PyFloat_FromDouble(_incs[6]),
		PyFloat_FromDouble(_incs[7]),
		PyFloat_FromDouble(_incs[8]),
		PyFloat_FromDouble(_incs[9]),
		PyFloat_FromDouble(_incs[10]),
		PyFloat_FromDouble(_incs[11]),
		PyFloat_FromDouble(_fval[0]),
		PyFloat_FromDouble(_fval[1]),
		PyFloat_FromDouble(_fval[2]),
		PyFloat_FromDouble(_fval[3])
	);

	PyObject* _result = PyObject_CallObject(_function, _args);
	result = PyFloat_AsDouble(_result);

	Py_Finalize();
	return result;
}

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