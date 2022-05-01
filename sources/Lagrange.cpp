#include "Lagrange.h"

Lagrange::Lagrange(const uint& dimension,
	const uint& constraints,
	Parameters& parameters,
	Problem& problem) :
	ParabolaMethod(dimension, constraints, parameters, problem) {}

void Lagrange::calculate_localLipshConst(const uint& id_hyp) {
	Hyperinterval& hyp1 = _intervals[id_hyp];
	Hyperinterval& hyp2 = _intervals[_generated_intervals - 2];
	Hyperinterval& hyp3 = _intervals[_generated_intervals - 1];
	
	double fa, fv, fu, fb;
	double dist = 0.0;
	for (uint i = 0; i < _dimension; ++i)
		dist = dist + ((double)_coords[hyp1.get_coordA() + i] - 
				       (double)_coords[hyp3.get_coordB() + i]) * 
				      ((double)_coords[hyp1.get_coordA() + i] - 
				       (double)_coords[hyp3.get_coordB() + i]);
	
	for (uint i = 0; i < _constraints + 1; ++i) {
		fa = _evaluations[hyp1.get_evalA() + i];
		fv = _evaluations[hyp2.get_evalA() + i];
		fu = _evaluations[hyp2.get_evalB() + i];
		fb = _evaluations[hyp3.get_evalB() + i];
		_localLipshEval[i] = 4 * (fb - fu - fv + fa) / dist;
		double conv = _localLipshEval[i] * (double)MAX_POWER_THREE * (double)MAX_POWER_THREE;
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