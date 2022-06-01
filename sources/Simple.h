#ifndef SIMPLE_H
#define SIMPLE_H

#include "Piyavskii.h"

class Simple : public Piyavskii {
private:
	void calculate_characteristic(const uint& id_hyp) override;
	void give_borders(double& l, double& r, Hyperinterval& hyp);
	void update_minimum(const FunctionsValues&, const uint&) override;
public:
	Simple(const uint& dimension,
		const uint& constraints,
		Parameters& parameters,
		Problem& problem);
	~Simple() {}
};

#endif // SIMPLE_H