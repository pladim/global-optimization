#ifndef DIVIDE_BY_THREE_H
#define DIVIDE_BY_THREE_H

#include <deque>

#include "synonymous_types.h"
#include "Problem.h"
#include "Hyperinterval.h"
#include "Point.h"
#include "Parameters.h"

enum class Mode { stop_by_precision, test };

class DivideByThree {
protected:
    // параметры задачи
    // размерность задачи
    uint _dimension;
    // число ограничений в задаче
    uint _constraints;
    // структура содержащая прочие параметры
    Parameters& _parameters;
    // декоратор для целевой функции и ограничений
    Problem& _problem;
    // дек координат порождённых точек
    std::deque<uint> _coords;
    // дек измерений в порождённых точках
    std::deque<FunctionValue> _evaluations;
    // дек порождённых гиперинтервалов
    std::deque<Hyperinterval> _intervals;
    // дек порождённых точек
    std::deque<Point> _points;
    // число сгенерированных точек и интервалов методом
    uint _generated_points;
    uint _generated_intervals;
    // для нужд метода
    EncodedCoordinates _transit1;
    EncodedCoordinates _transit2;
    // информация о текущем минимуме функции
    FunctionValue _current_minimum;
    uint _id_minimum;
    // какую ось будем делить
    uint _divided_axis;
    // максимальное число итераций
    uint _max_it;
    // число итераций
    uint _iteration;
    bool _solved;
    // вектор расстояний между оценками усл. глоб.-опт. значения с шага на шаг
    std::vector<double> _distances;
    // режим, при котором запускается метод
    Mode _state;
    // счётчик, отслеживающий сколько раз не менялся глобальный минимум
    uint _changes;
protected:
	// создать первый гиперинтервал
	void initialization();
	// поделить гиперинтервал на три части
	void trisect_interval(const uint& id_hyp);
	// заполнить гиперинтервал данными
	void fill_intervals(Hyperinterval& parent, const uint& id_u, const uint& id_v);
	// вычислить характеристику гиперинтервала
	virtual void calculate_characteristic(const uint& id_hyp) = 0;
    // выдать идентификатор точке
    uint generate_id();
    // выдать идентификатор гиперинтервалу
    uint generate_hyp();
    // вычислить значение функции в точке
    void compute_evaluations(const uint& id_point);
    // обновить оценку минимума
    virtual void update_minimum(const FunctionsValues&, const uint&);
    // вычислить длину диагонали гиперинтервала
    void compute_diagonal(const uint& id_hyp);
    // выбрать наилучший гиперинтервал для деления
    virtual uint optimal_to_trisect() = 0;
    // сделать шаг метода
    virtual uint iterate(const uint& id_hyp) = 0;
    // вычислить расстояние
    void calc_distance(const uint& idx);
protected:
    // методы расширения деков
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
    bool solved() const { return _solved; }
    std::vector<double> get_distances() { return _distances; }
    CoordinatesValues get_min_point();
    virtual ~DivideByThree() {}
};

#endif // DIVIDE_BY_THREE_H
