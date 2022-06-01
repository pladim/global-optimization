#include "Piyavskii.h"

Piyavskii::Piyavskii(const uint& dimension,
					 const uint& constraints,
					 Parameters& parameters,
					 Problem& problem) : 
	ParabolaMethod(dimension, constraints, parameters, problem) {
}

void Piyavskii::calculate_localLipshConst(const uint& id_hyp) {
	Hyperinterval& hyp = _intervals[id_hyp];

	for (uint i = 0; i < _constraints + 1; ++i) {
		_localLipshEval[i] = std::abs(_evaluations[hyp.get_evalB() + i] - _evaluations[hyp.get_evalA() + i]);
		_localLipshEval[i] /= hyp.get_diagonal();
	}

	hyp.update_localLipQueues(
		_localLipshEval,
		_parameters._delta / ((double)MAX_POWER_THREE * (double)MAX_POWER_THREE)
	);
}