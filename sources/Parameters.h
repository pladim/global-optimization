#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "synonymous_types.h"

struct Parameters {
	uint _dimension;
	uint _constraints;
	uint _queueDepth;
	uint _max_it;
	GainLipshConstant _gainLocalObj;
	GainLipshConstant _gainLocalCst;
	GainLipshConstant _gainGlobalObj;
	GainLipshConstant _gainGlobalCst;
	double _delta;
	double _criticalSize;
	double _eps;

	GainLipshConstant _Gain;
	GainLipshConstant _Reduce;

	double _Delta;
	uint _iter_thr;
};

#endif // PARAMETERS_H