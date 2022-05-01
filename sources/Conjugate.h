#ifndef CONJUGATE_H
#define CONJUGATE_H

#include "ParabolaMethod.h"

class Conjugate : public ParabolaMethod {
protected:
	void calculate_localLipshConst(const uint& id_hyp) override;
public:
	Conjugate(const uint& dimension,
			  const uint& constraints,
			  Parameters& parameters,
			  Problem& problem);
	virtual ~Conjugate() {}
};

#endif // CONJUGATE_H
