#ifndef TRANSFORMPMWITHSM_H
#define TRANSFORMPMWITHSM_H

#include "WithSimplexMethod.h"

class TransformPMwithSM : public WithSimplexMethod {
private:
	void calculate_characteristic(const uint& id_hyp) override;
	double golden_ratio(double a,
						double b,
						const double& e,
						const uint& id_hyp);
	double calculate_residual(const double& t, const double& e,
							  const uint& id_hyp);
	void update_minimum(const FunctionsValues&, const uint&) override;
public:
	TransformPMwithSM(const uint& dimension,
					  const uint& constraints,
					  Parameters& parameters,
					  Problem& problem);
	~TransformPMwithSM() {}
};

#endif // TRANSFORMPMWITHSM_H