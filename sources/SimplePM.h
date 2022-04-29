#ifndef SIMPLEPM_H
#define SIMPLEPM_H

#include "DivideByThree.h"

class SimplePM : public DivideByThree {
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
protected:
	void give_borders(double& l, double& r, Hyperinterval& hyp);
	void balance(double& _lipshConst) const;
private:
	void decode_and_save(const uint& pos, const uint& order);
	void projection(const uint& order,
					std::vector<double>& e1,
					std::vector<double>& e2);
	void calculate_and_project(const uint& axis);
	double scalar_product(const std::vector<double>::iterator& a,
						  const std::vector<double>::iterator& b);
	void generate_right_part(const uint& function,
							 const uint& eval_a,
							 const uint& eval_v,
							 const uint& eval_u,
							 const uint& eval_b);
	double get_solution();
public:
	SimplePM(const uint& dimension,
			 const uint& constraints,
			 Parameters& parameters,
			 Problem& problem);
	void solve() override;
	~SimplePM() {}
};

#endif SIMPLEPM_H
