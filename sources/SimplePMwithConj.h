#ifndef SIMPLEPMWITHCONJ_H
#define SIMPLEPMWITHCONJ_H

#include "Conjugate.h"

class SimplePMwithConj : public Conjugate {
private:
	void calculate_characteristic(const uint& id_hyp) override;
	void give_borders(double& l, double& r, Hyperinterval& hyp);
public:
	SimplePMwithConj(const uint& dimension,
					 const uint& constraints,
					 Parameters& parameters,
					 Problem& problem);
	~SimplePMwithConj() {}
};

#endif // SIMPLEPMWITHCONJ_H
