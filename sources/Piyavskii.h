#ifndef PIYAVSKII_H
#define PIYAVSKII_H

#include "ParabolaMethod.h"

class Piyavskii : public ParabolaMethod {
protected:
	void calculate_localLipshConst(const uint& id_hyp) override;
	uint iterate(const uint& id_hyp) override;
public:
	Piyavskii(const uint& dimension,
			  const uint& constraints,
			  Parameters& parameters,
			  Problem& problem);
	virtual ~Piyavskii() {}
};

#endif // PIYAVSKII_H
