#ifndef LAGRANGE_H
#define LAGRANGE_H

#include "ParabolaMethod.h"

class Lagrange : public ParabolaMethod {
protected:
	void calculate_localLipshConst(const uint& id_hyp) override;
public:
	Lagrange(const uint& dimension,
			  const uint& constraints,
			  Parameters& parameters,
			  Problem& problem);
	virtual ~Lagrange() {}
};

#endif // LAGRANGE_H
