/**
 * \file parametrized_line.hpp
 * \brief This file declares the class for representing parametrization
 *        of a line segment in 2D. This class inherits from the
 *        Abstract base class representing parametrized curves
 * @see abstract_parametrized_curve.hpp
 *
 *  This File is a part of the 2D-Parametric BEM package
 */

#ifndef PARAMETRIZEDLINEHPP
#define PARAMETRIZEDLINEHPP

#include <utility>

#include <Eigen/Dense>

#include "abstract_parametrized_curve.hpp"

class ParametrizedLine : public AbstractParametrizedCurve {
  public:
    using Point = Eigen::Vector2d;

    /**
     * Constructor with start and end points of the 2D line
     *
     * @param start Starting point for the line segment
     * @param end Ending point for the line segment
     */
    ParametrizedLine(Point, Point);

    /**
     * @see AbstractParametrizedCurve::ParameterRange()
     */
    std::pair<double,double> ParameterRange(void) const;

    /**
     * @see AbstractParametrizedCurve::operator() (double)
     */
    Eigen::Vector2d operator() (double) const;

    /**
     * @see AbstractParametrizedCurve::Derivative(double)
     */
    Eigen::Vector2d Derivative(double) const;

  private:
    /**
     * private const fields for storing the starting and
     * ending point of the line.
     */
    const Point start;
    const Point end;
};

#endif //PARAMETRIZEDLINEHPP
