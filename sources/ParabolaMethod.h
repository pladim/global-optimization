#ifndef PARABOLAMETHOD_H
#define PARABOLAMETHOD_H

#include "DivideByThree.h"

class ParabolaMethod : public DivideByThree {
protected:
	bool _areAllCharInfty;
	bool _doesGlobalChange;
	FunctionsValues _localLipshEval;
	FunctionsValues _globalLipshEval;
protected:
	virtual void calculate_localLipshConst(const uint& id_hyp) = 0;
	void calculate_globalLipshConst(const uint& id_hyp);
	double mixedLipEval(const Hyperinterval& hyp, const uint& i);
	void update_all_charact();
	uint optimal_to_trisect() override;
	uint iterate(const uint& id_hyp) override;
	void balance(double& _lipshConst) const;
	uint min_by_charact();
	uint max_by_length();
public:
	ParabolaMethod(const uint& dimension,
				   const uint& constraints,
				   Parameters& parameters,
				   Problem& problem);
	void solve() override;
	virtual ~ParabolaMethod() {}
};

#endif // PARABOLAMETHOD_H
