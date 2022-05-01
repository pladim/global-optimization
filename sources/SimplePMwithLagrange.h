#ifndef SIMPLEPMWITHLAGRANGE_H
#define SIMPLEPMWITHLAGRANGE_H

#include "Lagrange.h"

class SimplePMwithLagrange : public Lagrange {
private:
	void calculate_characteristic(const uint& id_hyp) override;
	void give_borders(double& l, double& r, Hyperinterval& hyp);
public:
	SimplePMwithLagrange(const uint& dimension,
						 const uint& constraints,
						 Parameters& parameters,
						 Problem& problem);
	~SimplePMwithLagrange() {}
};

#endif // SIMPLEPMWITHLAGRANGE_H
