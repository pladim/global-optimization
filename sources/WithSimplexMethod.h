#ifndef WITHSIMPLEXMETHOD_H
#define WITHSIMPLEXMETHOD_H

#include "ortools/linear_solver/linear_solver.h"
#include "ParabolaMethod.h"

class WithSimplexMethod : public ParabolaMethod {
protected:
	CoordinatesValues _point;
	std::vector<double> _incs;
	std::vector<double> _fval;
	std::vector<double> _non_proj_incs;
protected:
	void calculate_localLipshConst(const uint& id_hyp) override;
	double get_solution();
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
public:
	WithSimplexMethod(const uint& dimension,
					 const uint& constraints,
					 Parameters& parameters,
					 Problem& problem);
	virtual ~WithSimplexMethod() {}
};

#endif // WITHSIMPLEXMETHOD_H