#ifndef DIVIDE_BY_THREE_H
#define DIVIDE_BY_THREE_H

#include <deque>

#include "synonymous_types.h"
#include "Problem.h"
#include "Hyperinterval.h"
#include "Point.h"
#include "Parameters.h"

class DivideByThree {
protected:
    // ��������� ������
    // ����������� ������
    uint _dimension;
    // ����� ����������� � ������
    uint _constraints;
    // ��������� ���������� ������ ���������
    Parameters& _parameters;
    // ��������� ��� ������� ������� � �����������
    Problem& _problem;
    // ��� ��������� ���������� �����
    std::deque<uint> _coords;
    // ��� ��������� � ���������� ������
    std::deque<FunctionValue> _evaluations;
    // ��� ���������� ���������������
    std::deque<Hyperinterval> _intervals;
    // ��� ���������� �����
    std::deque<Point> _points;
    // ����� ��������������� ����� � ���������� �������
    uint _generated_points;
    uint _generated_intervals;
    // ��� ���� ������
    EncodedCoordinates _transit1;
    EncodedCoordinates _transit2;
    // ���������� � ������� �������� �������
    FunctionValue _current_minimum;
    uint _id_minimum;
    // ����� ��� ����� ������
    uint _divided_axis;
    // ������������ ����� ��������
    uint _max_it;
    // ����� ��������
    uint _iteration;
protected:
	// ������� ������ �������������
	void initialization();
	// �������� ������������� �� ��� �����
	void trisect_interval(const uint& id_hyp);
	// ��������� ������������� �������
	void fill_intervals(Hyperinterval& parent, const uint& id_u, const uint& id_v);
	// ��������� �������������� ��������������
	virtual void calculate_characteristic(const uint& id_hyp) = 0;
    // ������ ������������� �����
    uint generate_id();
    // ������ ������������� ��������������
    uint generate_hyp();
    // ��������� �������� ������� � �����
    void compute_evaluations(const uint& id_point);
    // �������� ������ ��������
    virtual void update_minimum(const FunctionsValues&, const uint&);
    // ��������� ����� ��������� ��������������
    void compute_diagonal(const uint& id_hyp);
    // ������� ��������� ������������� ��� �������
    virtual uint optimal_to_trisect() = 0;
    // ������� ��� ������
    virtual uint iterate(const uint& id_hyp) = 0;
protected:
    // ������ ���������� �����
    void resize_points_deque();
    void resize_coords_deque();
    void resize_evaluations_deque();
    void resize_intervals_deque();
public:
    DivideByThree(const uint& dimension,
                  const uint& constraints,
                  Parameters& parameters,
                  Problem& problem);
    virtual void solve() = 0;
    void write_generated_points();
    void write_generated_intervals();
    uint get_gen() const { return _generated_points; }
    double get_min() const { return _current_minimum; }
    virtual ~DivideByThree() {}
};

#endif // DIVIDE_BY_THREE_H
