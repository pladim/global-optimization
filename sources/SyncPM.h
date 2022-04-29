#ifndef SYNCPM_H
#define SYNCPM_H

#include <Python.h>
#include "DivideByThree.h"

class SyncPM : public DivideByThree {
protected:
	bool _areAllCharInfty;
	bool _doesGlobalChange;
	FunctionsValues _localLipshEval;
	FunctionsValues _globalLipshEval;
private:
	CoordinatesValues _point;
	std::vector<double> _incs;
	std::vector<double> _fval;
	std::vector<double> _non_proj_incs;
protected:
	virtual void calculate_localLipshConst(const uint& id_hyp);
	void calculate_globalLipshConst(const uint& id_hyp);
	double mixedLipEval(const Hyperinterval& hyp, const uint& i);
	void update_all_charact();
	uint optimal_to_trisect() override;
	uint iterate(const uint& id_hyp) override;
	void calculate_characteristic(const uint& id_hyp) override;
protected:
	void give_borders(double& l, double& r, Hyperinterval& hyp);
	void balance(double& _lipshConst) const;
private:
	double golden_ratio(double a, double b, const double& e, const uint& id_hyp);
	double calculate_residual(const double& t, const double& e, const uint& id_hyp);
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
	SyncPM(const uint& dimension,
		   const uint& constraints,
		   Parameters& parameters,
		   Problem& problem);
	void solve() override;
	~SyncPM() {}
};

#endif // SYNCPM_H
