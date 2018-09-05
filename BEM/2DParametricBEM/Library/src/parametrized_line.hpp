/*
 * This hpp file declares the class representing parametrization
 * of a line in 2-D which is inherited from the
 * Abstract base class representing parametrized curves
 * (See abstract_parametrized_curve.hpp)
 */

#ifndef PARAMETRIZEDLINEHPP
#define PARAMETRIZEDLINEHPP

#include <utility>

#include <Eigen/Dense>

#include "abstract_parametrized_curve.hpp"

class ParametrizedLine : public AbstractParametrizedCurve {
  public:
    using Point = Eigen::Vector2d;

    /* Constructor with start and end points of the 2D line*/
    ParametrizedLine(Point, Point);

    std::pair<double,double> ParameterRange(void) const;

    Point operator() (double t) const;

    Eigen::Vector2d Derivative(double t) const;
  private:
    /* private const fields for the starting and ending point of the line*/
    const Point start;
    const Point end;
};

#endif //PARAMETRIZEDLINEHPP
