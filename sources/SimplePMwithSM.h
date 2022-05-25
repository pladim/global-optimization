#ifndef SIMPLEPMWITHSM_H
#define SIMPLEPMWITHSM_H

#include "WithSimplexMethod.h"

class SimplePMwithSM : public WithSimplexMethod {
private:
	void calculate_characteristic(const uint& id_hyp) override;
	void give_borders(double& l, double& r, Hyperinterval& hyp);
	void update_minimum(const FunctionsValues&, const uint&) override;
public:
	SimplePMwithSM(const uint& dimension,
				   const uint& constraints,
				   Parameters& parameters,
				   Problem& problem);
	~SimplePMwithSM() {}
};

#endif // SIMPLEPMWITHSM_H
