#ifndef PROBLEM_H
#define PROBLEM_H

#include "synonymous_types.h"

class Problem {
	// размерность задачи и число ограничений
	uint _dimension;
	uint _constraints;
	// область поиска
	CoordinatesValues _left_borders;
	CoordinatesValues _right_borders;
	// для вычисления значений целевой функции и ограничений
	FunctionsCalculator _problem;
	// вектор для хранения расшифрованных координат
	CoordinatesValues _decoded_point;
	// 
	FunctionsValues _function_values;
	// истинная точка
	CoordinatesValues _true_min;
private:
	// вычислить функцию в точке для вещественных координат
	FunctionsValues& f(const CoordinatesValues& out);
public:
	Problem() {}
	Problem(const uint& dimension, 
			const uint& constraints,
			const CoordinatesValues& left_brd,
			const CoordinatesValues& right_brd,
			const FunctionsCalculator& problem);
	// расшифровать координаты
	CoordinatesValues& decode_coordinates(const EncodedCoordinates& out);
	CoordinatesValues& decode_coordinates01(const EncodedCoordinates& out);
	void set_true_minimum(const CoordinatesValues& out);
	CoordinatesValues get_true_minimum() { return _true_min; }
	// вычислить значений целевой функции и ограничений
	FunctionsValues& operator()(const EncodedCoordinates& out);
	~Problem() {}
};

#endif // PROBLEM_H
