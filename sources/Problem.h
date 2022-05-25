#ifndef PROBLEM_H
#define PROBLEM_H

#include "synonymous_types.h"

class Problem {
	// ����������� ������ � ����� �����������
	uint _dimension;
	uint _constraints;
	// ������� ������
	CoordinatesValues _left_borders;
	CoordinatesValues _right_borders;
	// ��� ���������� �������� ������� ������� � �����������
	FunctionsCalculator _problem;
	// ������ ��� �������� �������������� ���������
	CoordinatesValues _decoded_point;
	// 
	FunctionsValues _function_values;
	// �������� �����
	CoordinatesValues _true_min;
private:
	// ��������� ������� � ����� ��� ������������ ���������
	FunctionsValues& f(const CoordinatesValues& out);
public:
	Problem() {}
	Problem(const uint& dimension, 
			const uint& constraints,
			const CoordinatesValues& left_brd,
			const CoordinatesValues& right_brd,
			const FunctionsCalculator& problem);
	// ������������ ����������
	CoordinatesValues& decode_coordinates(const EncodedCoordinates& out);
	CoordinatesValues& decode_coordinates01(const EncodedCoordinates& out);
	void set_true_minimum(const CoordinatesValues& out);
	CoordinatesValues get_true_minimum() { return _true_min; }
	// ��������� �������� ������� ������� � �����������
	FunctionsValues& operator()(const EncodedCoordinates& out);
	~Problem() {}
};

#endif // PROBLEM_H
