#ifndef HYPERINTERVAL_H
#define HYPERINTERVAL_H

#include "dicretization.h"
#include "synonymous_types.h"

class Hyperinterval {
	// ������������� ��������������
	uint _id;

	// ������������� ����� ������� ���������
	uint _idA;
	uint _idB;

	// ���������� ����������� ������� ��������������
	uint _divisions;

	// �������������� ��������������
	double _charact;
	// ����� �������� ��������� ��������������
	double _diagonal;
	// �������, ������� ����� ����� ������������ � ���������� ��������� ���������� ���������
	std::vector<double> _add_const;

	// ������� ��������� ������ �������� �������
	std::vector<std::queue<LipschitzConstantValue>> _localLipEvaluations;
	// ������������ �������� ��������� ������ �������� �������
	std::vector<LipschitzConstantValue> _maxLipEvaluations;

	// ����������� ������
	static uint _dimension;
	// ���������� �����������
	static uint _constraints;
	// ������� �������
	static uint _queueDepth;

public:
	Hyperinterval();
	Hyperinterval(const Hyperinterval& hyperinterval);
	Hyperinterval& operator=(const Hyperinterval& hyperinterval);

	// ���������������� ����������� ����
	static void init_static(const uint&, const uint&, const uint&);

	// �������� ������������� ��������������
	uint get_id() const;
	// �������� ������������� ����� ������� ���������
	uint get_idA() const;
	uint get_idB() const;
	// �������� ������ �������� ��������� ����� ������� ���������
	uint get_coordA() const;
	uint get_coordB() const;
	// �������� ������ �������� ��������� � ������ ������� ���������
	uint get_evalA() const;
	uint get_evalB() const;
	// �������� �������� �������������� ��������������
	double get_charact(const bool& = false) const;
	// ��������� ���������� ���������
	void count_add_const(const double& ggobj, const double& ggcst);
	// �������� ���������� ���������
	double get_add_const(const uint& function) const;
	// ���������� ������������� ��������������
	void set_id(const uint&);
	// ���������� ������������� ����� ������� ���������
	void set_idA(const uint&);
	void set_idB(const uint&);
	// ���������� �������� ������� ��������� ��������������
	void set_diagonal(const double&);
	// ���������� �������� �������������� ��������������
	void set_charact(const double&);

	// ��������� ����� ������� ��������������
	void increase_divisions();
	// �������� ����� �������
	uint get_divisions() const;
	// �������� ����� ������������ ���, �� ������� ����� ������ �������������
	uint get_axis() const;
	// �������� ����� ������������ ���, �� ������� ���������� �������������
	uint get_previous_axis() const;
	// �������� ����� ����� ������������ ��� ��� ������� ��������������
	uint get_shift() const;
	// �������� ����� ������� ���������
	double get_diagonal() const;

	// ���������� ������ �������� �������, ��������� ����� �������, 
	// ���������� ������������ ������
	void update_localLipQueues(std::vector<LipschitzConstantValue>&, 
							   const double&);

	// ����� ������������ ������ ��������� �������
	void find_maxLipEval();

	// �������� �������� ������������ ������ �������� �������
	const std::vector<LipschitzConstantValue>& get_maxLipshEvaluations() const;
};

#endif HYPERINTERVAL_H
