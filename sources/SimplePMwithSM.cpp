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
	double result{ 0.0 };
	double norm21{ _incs[0] * _incs[0] + _incs[1] * _incs[1] };
	double norm31{ _incs[2] * _incs[2] + _incs[3] * _incs[3] };
	double norm41{ _incs[4] * _incs[4] + _incs[5] * _incs[5] };
	double norm32{ _incs[6] * _incs[6] + _incs[7] * _incs[7] };
	double norm42{ _incs[8] * _incs[8] + _incs[9] * _incs[9] };
	double norm43{ _incs[10] * _incs[10] + _incs[11] * _incs[11] };
	
	using namespace operations_research;
	std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GLOP"));
	const double infinity = solver->infinity();

	// variables
	MPVariable* const L   = solver->MakeNumVar(0.0, infinity, "L");
	MPVariable* const ga1 = solver->MakeNumVar(-infinity, infinity, "ga1");
	MPVariable* const ga2 = solver->MakeNumVar(-infinity, infinity, "ga2");
	MPVariable* const gv1 = solver->MakeNumVar(-infinity, infinity, "gv1");
	MPVariable* const gv2 = solver->MakeNumVar(-infinity, infinity, "gv2");
	MPVariable* const gu1 = solver->MakeNumVar(-infinity, infinity, "gu1");
	MPVariable* const gu2 = solver->MakeNumVar(-infinity, infinity, "gu2");
	MPVariable* const gb1 = solver->MakeNumVar(-infinity, infinity, "gb1");
	MPVariable* const gb2 = solver->MakeNumVar(-infinity, infinity, "gb2");

	// constraints, burst 1
	MPConstraint* const b1 = solver->MakeRowConstraint(-infinity, 0.0);
	b1->SetCoefficient(L, -norm21);
	b1->SetCoefficient(ga1, -_incs[0]);
	b1->SetCoefficient(ga2, -_incs[1]);
	b1->SetCoefficient(gv1, _incs[0]);
	b1->SetCoefficient(gv2, _incs[1]);
	//b1->SetCoefficient(gu1, 0);
	//b1->SetCoefficient(gu2, 0);
	//b1->SetCoefficient(gb1, 0);
	//b1->SetCoefficient(gb2, 0);

	MPConstraint* const b2 = solver->MakeRowConstraint(-infinity, 0.0);
	b2->SetCoefficient(L, -norm21);
	b2->SetCoefficient(ga1, _incs[0]);
	b2->SetCoefficient(ga2, _incs[1]);
	b2->SetCoefficient(gv1, -_incs[0]);
	b2->SetCoefficient(gv2, -_incs[1]);
	//b2->SetCoefficient(gu1, 0);
	//b2->SetCoefficient(gu2, 0);
	//b2->SetCoefficient(gb1, 0);
	//b2->SetCoefficient(gb2, 0);

	MPConstraint* const b3 = solver->MakeRowConstraint(-infinity, _fval[0] - _fval[1]);
	b3->SetCoefficient(L, -0.5 * norm21);
	b3->SetCoefficient(ga1, -_incs[0]);
	b3->SetCoefficient(ga2, -_incs[1]);
	//b3->SetCoefficient(gv1, 0);
	//b3->SetCoefficient(gv2, 0);
	//b3->SetCoefficient(gu1, 0);
	//b3->SetCoefficient(gu2, 0);
	//b3->SetCoefficient(gb1, 0);
	//b3->SetCoefficient(gb2, 0);

	MPConstraint* const b4 = solver->MakeRowConstraint(-infinity, _fval[1] - _fval[0]);
	b4->SetCoefficient(L, -0.5 * norm21);
	b4->SetCoefficient(ga1, _incs[0]);
	b4->SetCoefficient(ga2, _incs[1]);
	//b4->SetCoefficient(gv1, 0);
	//b4->SetCoefficient(gv2, 0);
	//b4->SetCoefficient(gu1, 0);
	//b4->SetCoefficient(gu2, 0);
	//b4->SetCoefficient(gb1, 0);
	//b4->SetCoefficient(gb2, 0);

	MPConstraint* const b5 = solver->MakeRowConstraint(-infinity, _fval[0] - _fval[1]);
	b5->SetCoefficient(L, -0.5 * norm21);
	//b5->SetCoefficient(ga1, 0);
	//b5->SetCoefficient(ga2, 0);
	b5->SetCoefficient(gv1, -_incs[0]);
	b5->SetCoefficient(gv2, -_incs[1]);
	//b5->SetCoefficient(gu1, 0);
	//b5->SetCoefficient(gu2, 0);
	//b5->SetCoefficient(gb1, 0);
	//b5->SetCoefficient(gb2, 0);

	MPConstraint* const b6 = solver->MakeRowConstraint(-infinity, _fval[1] - _fval[0]);
	b6->SetCoefficient(L, -0.5 * norm21);
	//b6->SetCoefficient(ga1, 0);
	//b6->SetCoefficient(ga2, 0);
	b6->SetCoefficient(gv1, _incs[0]);
	b6->SetCoefficient(gv2, _incs[1]);
	//b6->SetCoefficient(gu1, 0);
	//b6->SetCoefficient(gu2, 0);
	//b6->SetCoefficient(gb1, 0);
	//b6->SetCoefficient(gb2, 0);

	// constraints, burst 2
	MPConstraint* const b7 = solver->MakeRowConstraint(-infinity, 0.0);
	b7->SetCoefficient(L, -norm31);
	b7->SetCoefficient(ga1, -_incs[2]);
	b7->SetCoefficient(ga2, -_incs[3]);
	//b7->SetCoefficient(gv1, 0);
	//b7->SetCoefficient(gv2, 0);
	b7->SetCoefficient(gu1, _incs[2]);
	b7->SetCoefficient(gu2, _incs[3]);
	//b7->SetCoefficient(gb1, 0);
	//b7->SetCoefficient(gb2, 0);

	MPConstraint* const b8 = solver->MakeRowConstraint(-infinity, 0.0);
	b8->SetCoefficient(L, -norm31);
	b8->SetCoefficient(ga1, _incs[2]);
	b8->SetCoefficient(ga2, _incs[3]);
	//b8->SetCoefficient(gv1, 0);
	//b8->SetCoefficient(gv2, 0);
	b8->SetCoefficient(gu1, -_incs[2]);
	b8->SetCoefficient(gu2, -_incs[3]);
	//b8->SetCoefficient(gb1, 0);
	//b8->SetCoefficient(gb2, 0);

	MPConstraint* const b9 = solver->MakeRowConstraint(-infinity, _fval[0] - _fval[2]);
	b9->SetCoefficient(L, -0.5 * norm31);
	b9->SetCoefficient(ga1, -_incs[2]);
	b9->SetCoefficient(ga2, -_incs[3]);
	//b9->SetCoefficient(gv1, 0);
	//b9->SetCoefficient(gv2, 0);
	//b9->SetCoefficient(gu1, 0);
	//b9->SetCoefficient(gu2, 0);
	//b9->SetCoefficient(gb1, 0);
	//b9->SetCoefficient(gb2, 0);
	
	MPConstraint* const b10 = solver->MakeRowConstraint(-infinity, _fval[2] - _fval[0]);
	b10->SetCoefficient(L, -0.5 * norm31);
	b10->SetCoefficient(ga1, _incs[2]);
	b10->SetCoefficient(ga2, _incs[3]);
	//b10->SetCoefficient(gv1, 0);
	//b10->SetCoefficient(gv2, 0);
	//b10->SetCoefficient(gu1, 0);
	//b10->SetCoefficient(gu2, 0);
	//b10->SetCoefficient(gb1, 0);
	//b10->SetCoefficient(gb2, 0);

	MPConstraint* const b11 = solver->MakeRowConstraint(-infinity, _fval[0] - _fval[2]);
	b11->SetCoefficient(L, -0.5 * norm31);
	//b11->SetCoefficient(ga1, 0);
	//b11->SetCoefficient(ga2, 0);
	//b11->SetCoefficient(gv1, 0);
	//b11->SetCoefficient(gv2, 0);
	b11->SetCoefficient(gu1, -_incs[2]);
	b11->SetCoefficient(gu2, -_incs[3]);
	//b11->SetCoefficient(gb1, 0);
	//b11->SetCoefficient(gb2, 0);

	MPConstraint* const b12 = solver->MakeRowConstraint(-infinity, _fval[2] - _fval[0]);
	b12->SetCoefficient(L, -0.5 * norm31);
	//b12->SetCoefficient(ga1, 0);
	//b12->SetCoefficient(ga2, 0);
	//b12->SetCoefficient(gv1, 0);
	//b12->SetCoefficient(gv2, 0);
	b12->SetCoefficient(gu1, _incs[2]);
	b12->SetCoefficient(gu2, _incs[3]);
	//b12->SetCoefficient(gb1, 0);
	//b12->SetCoefficient(gb2, 0);

	// constraints, burst 3
	MPConstraint* const b13 = solver->MakeRowConstraint(-infinity, 0.0);
	b13->SetCoefficient(L, -norm41);
	b13->SetCoefficient(ga1, -_incs[4]);
	b13->SetCoefficient(ga2, -_incs[5]);
	//b13->SetCoefficient(gv1, 0);
	//b13->SetCoefficient(gv2, 0);
	//b13->SetCoefficient(gu1, 0);
	//b13->SetCoefficient(gu2, 0);
	b13->SetCoefficient(gb1, _incs[4]);
	b13->SetCoefficient(gb2, _incs[5]);

	MPConstraint* const b14 = solver->MakeRowConstraint(-infinity, 0.0);
	b14->SetCoefficient(L, -norm41);
	b14->SetCoefficient(ga1, _incs[4]);
	b14->SetCoefficient(ga2, _incs[5]);
	//b14->SetCoefficient(gv1, 0);
	//b14->SetCoefficient(gv2, 0);
	//b14->SetCoefficient(gu1, 0);
	//b14->SetCoefficient(gu2, 0);
	b14->SetCoefficient(gb1, -_incs[4]);
	b14->SetCoefficient(gb2, -_incs[5]);

	MPConstraint* const b15 = solver->MakeRowConstraint(-infinity, _fval[0] - _fval[4]);
	b15->SetCoefficient(L, -0.5 * norm41);
	b15->SetCoefficient(ga1, -_incs[4]);
	b15->SetCoefficient(ga2, -_incs[5]);
	//b15->SetCoefficient(gv1, 0);
	//b15->SetCoefficient(gv2, 0);
	//b15->SetCoefficient(gu1, 0);
	//b15->SetCoefficient(gu2, 0);
	//b15->SetCoefficient(gb1, 0);
	//b15->SetCoefficient(gb2, 0);

	MPConstraint* const b16 = solver->MakeRowConstraint(-infinity, _fval[4] - _fval[0]);
	b16->SetCoefficient(L, -0.5 * norm41);
	b16->SetCoefficient(ga1, _incs[4]);
	b16->SetCoefficient(ga2, _incs[5]);
	//b16->SetCoefficient(gv1, 0);
	//b16->SetCoefficient(gv2, 0);
	//b16->SetCoefficient(gu1, 0);
	//b16->SetCoefficient(gu2, 0);
	//b16->SetCoefficient(gb1, 0);
	//b16->SetCoefficient(gb2, 0);

	MPConstraint* const b17 = solver->MakeRowConstraint(-infinity, _fval[0] - _fval[4]);
	b17->SetCoefficient(L, -0.5 * norm41);
	//b17->SetCoefficient(ga1, 0);
	//b17->SetCoefficient(ga2, 0);
	//b17->SetCoefficient(gv1, 0);
	//b17->SetCoefficient(gv2, 0);
	//b17->SetCoefficient(gu1, 0);
	//b17->SetCoefficient(gu2, 0);
	b17->SetCoefficient(gb1, -_incs[4]);
	b17->SetCoefficient(gb2, -_incs[5]);

	MPConstraint* const b18 = solver->MakeRowConstraint(-infinity, _fval[4] - _fval[0]);
	b18->SetCoefficient(L, -0.5 * norm41);
	//b18->SetCoefficient(ga1, 0);
	//b18->SetCoefficient(ga2, 0);
	//b18->SetCoefficient(gv1, 0);
	//b18->SetCoefficient(gv2, 0);
	//b18->SetCoefficient(gu1, 0);
	//b18->SetCoefficient(gu2, 0);
	b18->SetCoefficient(gb1, _incs[4]);
	b18->SetCoefficient(gb2, _incs[5]);

	// constraints, burst 4
	MPConstraint* const b19 = solver->MakeRowConstraint(-infinity, 0.0);
	b19->SetCoefficient(L, -norm32);
	//b19->SetCoefficient(ga1, 0);
	//b19->SetCoefficient(ga2, 0);
	b19->SetCoefficient(gv1, -_incs[6]);
	b19->SetCoefficient(gv2, -_incs[7]);
	b19->SetCoefficient(gu1, _incs[6]);
	b19->SetCoefficient(gu2, _incs[7]);
	//b19->SetCoefficient(gb1, 0);
	//b19->SetCoefficient(gb2, 0);

	MPConstraint* const b20 = solver->MakeRowConstraint(-infinity, 0.0);
	b20->SetCoefficient(L, -norm32);
	//b20->SetCoefficient(ga1, 0);
	//b20->SetCoefficient(ga2, 0);
	b20->SetCoefficient(gv1, _incs[6]);
	b20->SetCoefficient(gv2, _incs[7]);
	b20->SetCoefficient(gu1, -_incs[6]);
	b20->SetCoefficient(gu2, -_incs[7]);
	//b20->SetCoefficient(gb1, 0);
	//b20->SetCoefficient(gb2, 0);

	MPConstraint* const b21 = solver->MakeRowConstraint(-infinity, _fval[1] - _fval[2]);
	b21->SetCoefficient(L, -0.5 * norm32);
	//b21->SetCoefficient(ga1, 0);
	//b21->SetCoefficient(ga2, 0);
	b21->SetCoefficient(gv1, -_incs[6]);
	b21->SetCoefficient(gv2, -_incs[7]);
	//b21->SetCoefficient(gu1, 0);
	//b21->SetCoefficient(gu2, 0);
	//b21->SetCoefficient(gb1, 0);
	//b21->SetCoefficient(gb2, 0);

	MPConstraint* const b22 = solver->MakeRowConstraint(-infinity, _fval[2] - _fval[1]);
	b22->SetCoefficient(L, -0.5 * norm32);
	//b22->SetCoefficient(ga1, 0);
	//b22->SetCoefficient(ga2, 0);
	b22->SetCoefficient(gv1, _incs[6]);
	b22->SetCoefficient(gv2, _incs[7]);
	//b22->SetCoefficient(gu1, 0);
	//b22->SetCoefficient(gu2, 0);
	//b22->SetCoefficient(gb1, 0);
	//b22->SetCoefficient(gb2, 0);

	MPConstraint* const b23 = solver->MakeRowConstraint(-infinity, _fval[1] - _fval[2]);
	b23->SetCoefficient(L, -0.5 * norm32);
	//b23->SetCoefficient(ga1, 0);
	//b23->SetCoefficient(ga2, 0);
	//b23->SetCoefficient(gv1, 0);
	//b23->SetCoefficient(gv2, 0);
	b23->SetCoefficient(gu1, -_incs[6]);
	b23->SetCoefficient(gu2, -_incs[7]);
	//b23->SetCoefficient(gb1, 0);
	//b23->SetCoefficient(gb2, 0);

	MPConstraint* const b24 = solver->MakeRowConstraint(-infinity, _fval[2] - _fval[1]);
	b24->SetCoefficient(L, -0.5 * norm32);
	//b24->SetCoefficient(ga1, 0);
	//b24->SetCoefficient(ga2, 0);
	//b24->SetCoefficient(gv1, 0);
	//b24->SetCoefficient(gv2, 0);
	b24->SetCoefficient(gu1, _incs[6]);
	b24->SetCoefficient(gu2, _incs[7]);
	//b24->SetCoefficient(gb1, 0);
	//b24->SetCoefficient(gb2, 0);

	// constraints, burst 5
	MPConstraint* const b25 = solver->MakeRowConstraint(-infinity, 0.0);
	b25->SetCoefficient(L, -norm42);
	//b25->SetCoefficient(ga1, 0);
	//b25->SetCoefficient(ga2, 0);
	b25->SetCoefficient(gv1, -_incs[8]);
	b25->SetCoefficient(gv2, -_incs[9]);
	//b25->SetCoefficient(gu1, 0);
	//b25->SetCoefficient(gu2, 0);
	b25->SetCoefficient(gb1, _incs[8]);
	b25->SetCoefficient(gb2, _incs[9]);

	MPConstraint* const b26 = solver->MakeRowConstraint(-infinity, 0.0);
	b26->SetCoefficient(L, -norm42);
	//b26->SetCoefficient(ga1, 0);
	//b26->SetCoefficient(ga2, 0);
	b26->SetCoefficient(gv1, _incs[8]);
	b26->SetCoefficient(gv2, _incs[9]);
	//b26->SetCoefficient(gu1, 0);
	//b26->SetCoefficient(gu2, 0);
	b26->SetCoefficient(gb1, -_incs[8]);
	b26->SetCoefficient(gb2, -_incs[9]);

	MPConstraint* const b27 = solver->MakeRowConstraint(-infinity, _fval[1] - _fval[3]);
	b27->SetCoefficient(L, -0.5 * norm42);
	//b27->SetCoefficient(ga1, 0);
	//b27->SetCoefficient(ga2, 0);
	b27->SetCoefficient(gv1, -_incs[8]);
	b27->SetCoefficient(gv2, -_incs[9]);
	//b27->SetCoefficient(gu1, 0);
	//b27->SetCoefficient(gu2, 0);
	//b27->SetCoefficient(gb1, 0);
	//b27->SetCoefficient(gb2, 0);

	MPConstraint* const b28 = solver->MakeRowConstraint(-infinity, _fval[3] - _fval[1]);
	b28->SetCoefficient(L, -0.5 * norm42);
	//b28->SetCoefficient(ga1, 0);
	//b28->SetCoefficient(ga2, 0);
	b28->SetCoefficient(gv1, _incs[8]);
	b28->SetCoefficient(gv2, _incs[9]);
	//b28->SetCoefficient(gu1, 0);
	//b28->SetCoefficient(gu2, 0);
	//b28->SetCoefficient(gb1, 0);
	//b28->SetCoefficient(gb2, 0);

	MPConstraint* const b29 = solver->MakeRowConstraint(-infinity, _fval[1] - _fval[3]);
	b29->SetCoefficient(L, -0.5 * norm42);
	//b29->SetCoefficient(ga1, 0);
	//b29->SetCoefficient(ga2, 0);
	//b29->SetCoefficient(gv1, 0);
	//b29->SetCoefficient(gv2, 0);
	//b29->SetCoefficient(gu1, 0);
	//b29->SetCoefficient(gu2, 0);
	b29->SetCoefficient(gb1, -_incs[8]);
	b29->SetCoefficient(gb2, -_incs[9]);

	MPConstraint* const b30 = solver->MakeRowConstraint(-infinity, _fval[3] - _fval[1]);
	b30->SetCoefficient(L, -0.5 * norm42);
	//b30->SetCoefficient(ga1, 0);
	//b30->SetCoefficient(ga2, 0);
	//b30->SetCoefficient(gv1, 0);
	//b30->SetCoefficient(gv2, 0);
	//b30->SetCoefficient(gu1, 0);
	//b30->SetCoefficient(gu2, 0);
	b30->SetCoefficient(gb1, _incs[8]);
	b30->SetCoefficient(gb2, _incs[9]);

	// constraints, burst 6
	MPConstraint* const b31 = solver->MakeRowConstraint(-infinity, 0.0);
	b31->SetCoefficient(L, -norm43);
	//b31->SetCoefficient(ga1, 0);
	//b31->SetCoefficient(ga2, 0);
	//b31->SetCoefficient(gv1, 0);
	//b31->SetCoefficient(gv2, 0);
	b31->SetCoefficient(gu1, -_incs[10]);
	b31->SetCoefficient(gu2, -_incs[11]);
	b31->SetCoefficient(gb1, _incs[10]);
	b31->SetCoefficient(gb2, _incs[11]);

	MPConstraint* const b32 = solver->MakeRowConstraint(-infinity, 0.0);
	b32->SetCoefficient(L, -norm43);
	//b32->SetCoefficient(ga1, 0);
	//b32->SetCoefficient(ga2, 0);
	//b32->SetCoefficient(gv1, 0);
	//b32->SetCoefficient(gv2, 0);
	b32->SetCoefficient(gu1, _incs[10]);
	b32->SetCoefficient(gu2, _incs[11]);
	b32->SetCoefficient(gb1, -_incs[10]);
	b32->SetCoefficient(gb2, -_incs[11]);

	MPConstraint* const b33 = solver->MakeRowConstraint(-infinity, _fval[2] - _fval[3]);
	b33->SetCoefficient(L, -0.5 * norm43);
	//b33->SetCoefficient(ga1, 0);
	//b33->SetCoefficient(ga2, 0);
	//b33->SetCoefficient(gv1, 0);
	//b33->SetCoefficient(gv2, 0);
	b33->SetCoefficient(gu1, -_incs[10]);
	b33->SetCoefficient(gu2, -_incs[11]);
	//b33->SetCoefficient(gb1, 0);
	//b33->SetCoefficient(gb2, 0);

	MPConstraint* const b34 = solver->MakeRowConstraint(-infinity, _fval[3] - _fval[2]);
	b34->SetCoefficient(L, -0.5 * norm43);
	//b34->SetCoefficient(ga1, 0);
	//b34->SetCoefficient(ga2, 0);
	//b34->SetCoefficient(gv1, 0);
	//b34->SetCoefficient(gv2, 0);
	b34->SetCoefficient(gu1, _incs[10]);
	b34->SetCoefficient(gu2, _incs[11]);
	//b34->SetCoefficient(gb1, 0);
	//b34->SetCoefficient(gb2, 0);

	MPConstraint* const b35 = solver->MakeRowConstraint(-infinity, _fval[2] - _fval[3]);
	b35->SetCoefficient(L, -0.5 * norm43);
	//b35->SetCoefficient(ga1, 0);
	//b35->SetCoefficient(ga2, 0);
	//b35->SetCoefficient(gv1, 0);
	//b35->SetCoefficient(gv2, 0);
	//b35->SetCoefficient(gu1, 0);
	//b35->SetCoefficient(gu2, 0);
	b35->SetCoefficient(gb1, -_incs[10]);
	b35->SetCoefficient(gb2, -_incs[11]);

	MPConstraint* const b36 = solver->MakeRowConstraint(-infinity, _fval[3] - _fval[2]);
	b36->SetCoefficient(L, -0.5 * norm43);
	//b36->SetCoefficient(ga1, 0);
	//b36->SetCoefficient(ga2, 0);
	//b36->SetCoefficient(gv1, 0);
	//b36->SetCoefficient(gv2, 0);
	//b36->SetCoefficient(gu1, 0);
	//b36->SetCoefficient(gu2, 0);
	b36->SetCoefficient(gb1, _incs[10]);
	b36->SetCoefficient(gb2, _incs[11]);

	// objective function
	MPObjective* const objective = solver->MutableObjective();
	objective->SetCoefficient(L, 1);
	//objective->SetCoefficient(ga1, 0);
	//objective->SetCoefficient(ga2, 0);
	//objective->SetCoefficient(gv1, 0);
	//objective->SetCoefficient(gv2, 0);
	//objective->SetCoefficient(gu1, 0);
	//objective->SetCoefficient(gu2, 0);
	//objective->SetCoefficient(gb1, 0);
	//objective->SetCoefficient(gb2, 0);
	objective->SetMinimization();

	solver->Solve();
	result = objective->Value();

	/*Py_Initialize();
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

	Py_Finalize();*/
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