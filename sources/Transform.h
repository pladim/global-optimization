#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "Piyavskii.h"

class Transform : public Piyavskii {
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
	Transform(const uint& dimension,
			  const uint& constraints,
			  Parameters& parameters,
			  Problem& problem);
	~Transform() {}
};

#endif // TRANSFORM_H